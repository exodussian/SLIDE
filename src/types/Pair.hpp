/*
 * Pair.hpp
 *
 * This creates a unified input/output system for values
 * coming in pairs for both positive/negative electrodes.
 *  Created on: 19 Jun 2023
 *   Author(s): Volkan Kumtepeli
 */

namespace slide {

template <typename T = double>
struct ValuePair
{
  T pos{}, neg{}; //<! A value for positive/negative electrode
  [[nodiscard]] constexpr auto begin() noexcept { return &pos; }
  [[nodiscard]] constexpr auto end() noexcept { return &neg + 1; }

  auto &operator[](size_t i) { return *(begin() + i); }
  auto operator[](size_t i) const { return *(begin() + i); }
};

}; // namespace slide
