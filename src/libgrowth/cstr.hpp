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

#ifndef _CSTR_H
#define _CSTR_H

template <unsigned char... Chars> struct string_t
{
    static constexpr unsigned int size() { return sizeof...(Chars); }
    static const char *data()
    {
        static constexpr const char data[] = {Chars...};
        return data;
    }
};

namespace detail
{
template <typename Str, unsigned int N, unsigned char... Chars>
struct make_string_t : make_string_t<Str, N - 1, Str().chars[N - 1], Chars...>
{
};

template <typename Str, unsigned char... Chars>
struct make_string_t<Str, 0, Chars...>
{
    typedef string_t<Chars...> type;
};
} // namespace detail

#define CSTR(str)                                                              \
    [] {                                                                       \
        struct Str                                                             \
        {                                                                      \
            const char *chars = str;                                           \
        };                                                                     \
        return detail::make_string_t<Str, sizeof(str)>::type();                \
    }()

#endif
