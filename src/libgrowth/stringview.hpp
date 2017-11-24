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

#ifndef _STRINGVIEW_H
#define _STRINGVIEW_H

#include <cstring>

struct stringview {
  template <unsigned int N>
  constexpr stringview(const char (&ar)[N]) : begin(ar),size((ar[N-1]==0) ? N-1 : N) {}  // strips trailing \0     // implicit

  template <typename String,typename Sfinae=decltype(((String*)0)->c_str(),((String*)0)->size())>
  constexpr stringview(String&& str) : begin(str.c_str()),size(str.size()) {}

  stringview(const char *begin) : begin(begin),size(std::strlen(begin)) {}

  constexpr stringview(const char *begin,unsigned int size) : begin(begin),size(size) {}

  constexpr bool empty() const {
    return (size==0);
  }

  constexpr char operator*() const {
    // assert(!empty());  // or: throw ?
    return *begin;
  }

  constexpr stringview substr(unsigned int start) const {
    return { begin+start,
             (start<size) ? size-start : 0 };
  }

  constexpr stringview substr(unsigned int start,unsigned int len) const {
    return { begin+start,
             (start<size) ?
               (len<size-start) ? len : size-start
             : 0 };
  }

private:
  const char *begin;
  unsigned int size;
};

#endif
