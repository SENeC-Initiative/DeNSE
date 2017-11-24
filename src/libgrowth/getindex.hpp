/*
 * Copyright (c) 2016 Tobias Hoffmann
 *
 * License: http://opensource.org/licenses/MIT
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 */

#ifndef _GETINDEX_H
#define _GETINDEX_H

// provides  pack_tools::get_index<I>(Ts&&... ts)

// (similar to what  std::get<I>(std::make_tuple(ts...))  does)

namespace pack_tools {

namespace detail {
  template <unsigned int> struct int_c {};

  template <unsigned int I>
  constexpr void *get_index_impl(int_c<I>) // invalid index
  {
    return {};
  }

  template <typename T0,typename... Ts>
  constexpr T0&& get_index_impl(int_c<0>,T0&& t0,Ts&&... ts)
  {
    return (T0&&)t0;
  }

  template <unsigned int I,typename T0,typename... Ts>
  constexpr auto get_index_impl(int_c<I>,T0&& t0,Ts&&... ts)
    -> decltype(get_index_impl(int_c<I-1>(),(Ts&&)ts...))
  {
    return get_index_impl(int_c<I-1>(),(Ts&&)ts...);
  }
} // namespace detail

template <unsigned int I,typename... Ts>
constexpr auto get_index(Ts&&... ts)
  -> decltype(detail::get_index_impl(detail::int_c<I>(),(Ts&&)ts...))
{
  static_assert((I<sizeof...(Ts)),"Invalid Index");
  return detail::get_index_impl(detail::int_c<I>(),(Ts&&)ts...);
}

} // namespace pack_tools

#endif
