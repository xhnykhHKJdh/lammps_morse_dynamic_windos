/* ----------------------------------------------------------------------
   Custom Morse with dynamic contour window to enforce 1D sliding
   along a single DNA contour for one specific cohesin head.

   See header for usage.
------------------------------------------------------------------------- */

#include "pair_morse_dynamic_window.h"
#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "neigh_list.h"
#include "utils.h"
#include "neighbor.h"
#include "neigh_list.h"
#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
/* ---------------------------------------------------------------------- */

PairMorseDynamicWindow::PairMorseDynamicWindow(LAMMPS *lmp)
    : PairMorse(lmp)
{
  window_half_      = 0;
  right_head_tag_   = 0;
  prev_closest_tag_ = 0;
}

/* ---------------------------------------------------------------------- */
/*  pair_style morse/dynamic_window cut [window_half] [head_tag]
 *
 *  - cut         : 全局 cutoff，同普通 morse
 *  - window_half : 允许的 tag 范围半宽（整数，单位 = beads 个数），
 *                  例如 5 表示只允许 [tag* - 5, tag* + 5]
 *  - head_tag    : 需要施加 1D 滑动约束的 cohesin 头的 tag（比如 1202）
 *
 *  如果 window_half_ == 0 或 head_tag_ == 0，则退化为普通 Morse。
 * ---------------------------------------------------------------------- */

void PairMorseDynamicWindow::settings(int narg, char **arg)
{
  // 允许 1~3 个参数：cut, [window_half], [head_tag]
  if (narg < 1 || narg > 3)
    error->all(FLERR, "Illegal pair_style morse/dynamic_window command");

  // 第一个参数：全局 cutoff（完全照抄 PairMorse::settings 的写法）
  cut_global = utils::numeric(FLERR, arg[0], false, lmp);

  window_half_      = 0;
  right_head_tag_   = 0;
  prev_closest_tag_ = 0;

  // 第二个参数：窗口半宽（可以是浮点，内部四舍五入成 int）
  if (narg >= 2) {
    double w = utils::numeric(FLERR, arg[1], false, lmp);
    if (w < 0.0)
      error->all(FLERR,
                 "pair_style morse/dynamic_window: window_half must be >= 0");
    window_half_ = static_cast<int>(w + 0.5);
  }

  // 第三个参数：cohesin 头的 tag
  if (narg >= 3) {
    right_head_tag_ = utils::tnumeric(FLERR, arg[2], false, lmp);
  }

  // 像原始 PairMorse 一样，如果已经分配过，就重置所有已设置的 cut[i][j]
  if (allocated) {
    for (int i = 1; i <= atom->ntypes; i++)
      for (int j = i; j <= atom->ntypes; j++)
        if (setflag[i][j]) cut[i][j] = cut_global;
  }
}

/* ---------------------------------------------------------------------- */
/*  compute:
 *  和原始 PairMorse 几乎一样，只在「cohesin-head ↔ DNA」那一类 pair 上
 *  做一个 tag-based window 的筛选：
 *
 *  - 通过 right_head_tag_ 找到那个特定的 cohesin 头；
 *  - 只允许与 prev_closest_tag_ ± window_half_ 范围内的 DNA tag 发生 Morse；
 *  - 同时，在这一 timestep 内记录距离这个头最近的 DNA 的 tag，
 *    在循环结束后更新 prev_closest_tag_，实现沿 contour 的“自然滑动”。
 * ---------------------------------------------------------------------- */

void PairMorseDynamicWindow::compute(int eflag, int vflag)
{
  int i, j, ii, jj, inum, jnum, itype, jtype;
  double xtmp, ytmp, ztmp, delx, dely, delz, evdwl, fpair;
  double rsq, r, dr, dexp, factor_lj;
  int *ilist, *jlist, *numneigh, **firstneigh;

  evdwl = 0.0;
  ev_init(eflag, vflag);

  double **x    = atom->x;
  double **f    = atom->f;
  int *type     = atom->type;
  tagint *tag   = atom->tag;
  int nlocal    = atom->nlocal;
  double *special_lj = force->special_lj;
  int newton_pair    = force->newton_pair;

  inum       = list->inum;
  ilist      = list->ilist;
  numneigh   = list->numneigh;
  firstneigh = list->firstneigh;

  // 这一步中找到的「最近 DNA tag」
  tagint current_closest_tag = prev_closest_tag_;
  double min_rsquared        = 1.0e30;

  // 只有当我们已经有合法的 prev_closest_tag_ 时，窗口才真正生效
  const bool window_active =
      (window_half_ > 0 && right_head_tag_ > 0 && prev_closest_tag_ > 0);

  // --------------------------------------------------
  // 主循环：几乎照抄原始 PairMorse::compute
  // --------------------------------------------------

  for (ii = 0; ii < inum; ii++) {
    i    = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum  = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j         = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      j        &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq  = delx * delx + dely * dely + delz * delz;
      jtype = type[j];

      if (rsq >= cutsq[itype][jtype]) continue;

      // -----------------------------
      // 这里做“拓扑窗口”的筛选
      // -----------------------------
      bool allow_morse = true;

      if (right_head_tag_ > 0) {
        tagint ti = tag[i];
        tagint tj = tag[j];

        bool has_head =
            (ti == right_head_tag_) || (tj == right_head_tag_);

        if (has_head) {
          // 哪一个是 DNA 的 tag？（假设只有一个是 head，另一个是 DNA）
          tagint dna_tag = (ti == right_head_tag_) ? tj : ti;

          // 用真实空间距离来更新“最近 DNA”
          if (rsq < min_rsquared) {
            min_rsquared        = rsq;
            current_closest_tag = dna_tag;
          }

          // 如果窗口已经激活，就只允许窗口内的 bead
          if (window_active) {
            tagint lo = prev_closest_tag_ - window_half_;
            tagint hi = prev_closest_tag_ + window_half_;
            if (dna_tag < lo || dna_tag > hi) {
              allow_morse = false;  // 这里就是禁止“桥接”的关键一步
            }
          }
        }
      }

      if (!allow_morse) continue;

      // -----------------------------
      // 正常的 Morse 计算（照抄 PairMorse）
      // -----------------------------

      r    = sqrt(rsq);
      dr   = r - r0[itype][jtype];
      dexp = exp(-alpha[itype][jtype] * dr);
      fpair = factor_lj * morse1[itype][jtype] * (dexp * dexp - dexp) / r;

      f[i][0] += delx * fpair;
      f[i][1] += dely * fpair;
      f[i][2] += delz * fpair;
      if (newton_pair || j < nlocal) {
        f[j][0] -= delx * fpair;
        f[j][1] -= dely * fpair;
        f[j][2] -= delz * fpair;
      }

      if (eflag) {
        evdwl = d0[itype][jtype] * (dexp * dexp - 2.0 * dexp) -
                offset[itype][jtype];
        evdwl *= factor_lj;
      }

      if (evflag)
        ev_tally(i, j, nlocal, newton_pair, evdwl, 0.0, fpair,
                 delx, dely, delz);
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();

  // --------------------------------------------------
  // 循环结束后，更新 prev_closest_tag_
  // 第一个 timestep：prev_closest_tag_ = 0，所以窗口还没激活，
  // 这一步会用“最近的 DNA”初始化它，后续才开始真正施加窗口。
  // --------------------------------------------------
  if (current_closest_tag > 0 &&
      current_closest_tag != prev_closest_tag_) {
    prev_closest_tag_ = current_closest_tag;
  }
}
