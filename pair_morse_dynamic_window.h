/* ----------------------------------------------------------------------
   Custom Morse with dynamic contour window to enforce 1D sliding
   along a single DNA contour (no bridging to other segments).

   pair_style morse/dynamic_window cut [window_half] [head_tag]

   - cut         : global cutoff (same as normal morse)
   - window_half : half-width of allowed DNA tag window (integer, in beads)
                   if 0, behaves like normal Morse (no window)
   - head_tag    : atom tag of the sliding cohesin head (e.g. 1202)
                   if 0, behaves like normal Morse

   Example (hybrid/overlay with WCA background):

     pair_style  hybrid/overlay lj/cut 1.122462048 morse/dynamic_window 1.0 5 1202
     pair_coeff  1 1 lj/cut 1.0 1.0 1.122462048
     pair_coeff  1 3 morse/dynamic_window 6.0 30.0 0.0 1.0

------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// 注册新的 pair_style 名字：morse/dynamic_window
PairStyle(morse/dynamic_window, PairMorseDynamicWindow);
#else

#ifndef LMP_PAIR_MORSE_DYNAMIC_WINDOW_H
#define LMP_PAIR_MORSE_DYNAMIC_WINDOW_H

#include "pair_morse.h"

namespace LAMMPS_NS {

class PairMorseDynamicWindow : public PairMorse {
 public:
  PairMorseDynamicWindow(class LAMMPS *);
  ~PairMorseDynamicWindow() override = default;

  void compute(int, int) override;
  void settings(int, char **) override;

 protected:
  int    window_half_;       // contour 窗口半宽（bead 个数）
  tagint right_head_tag_;    // 要施加拓扑约束的 cohesin 头的 tag
  tagint prev_closest_tag_;  // 上一 timestep 最近的 DNA bead 的 tag
};

}  // namespace LAMMPS_NS

#endif
#endif
