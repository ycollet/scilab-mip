/**********************************************************************
 * tokenize() - Tokenizes a sequence into parts by separators.        *
 * Copyright (C) 2000  C. Brandon Forehand                            *
 * <b4hand@users.sourceforge.net>                                     *
 *                                                                    *
 * This code is free software; you can redistribute it and/or         *
 * modify it under the terms of the GNU General Public License        *
 * as published by the Free Software Foundation; either version 2     *
 * of the License, or (at your option) any later version.             *
 *                                                                    *
 * This program is distributed in the hope that it will be useful,    *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of     *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the      *
 * GNU General Public License for more details.                       *
 **********************************************************************/

#ifndef TOKENIZE_H
#define TOKENIZE_H

#include <algorithm>
#include <cctype>
#include <iterator>
#include <functional>
#include <string>
#include <vector>

/* KLUDGE: ANSI/ISO C but not ANSI C++ allows for certain functions to
   be shadowed by macros, but in GCC 2.95.x this causes a parse error
   in is_space.  This allows code to be compiled in both GCC 3.0 and
   GCC 2.95.x.  */
#ifdef isspace
#undef isspace
#endif

//------------------------------------------------------------------------------
/** This is the default predicate for the tokenize function.  It is necessary
    because the prototype for isspace() is int isspace(int). */
//------------------------------------------------------------------------------
#if !defined(WIN32)
struct is_space : public std::unary_function<char,bool>
{
    bool operator() (char c) const
        { return std::isspace(c) != 0; }
};
#else
struct is_space : public std::unary_function<char,bool>
{
    bool operator() (char c) const
        { return isspace(c) != 0; }
};
#endif

//------------------------------------------------------------------------------
/** Powerful general predicate for tokenize().  The predicate returns
    true for all arguments that are contained in the sequence passed
    in at construction.  For example, the predicate created by
    recognize(" \v\f\n\r\t") is the same as is_space() above. */
//------------------------------------------------------------------------------
template<class Container = std::string>
class recognize_t : 
    public std::unary_function<typename Container::value_type,bool>
{
public:
    recognize_t(const Container &c) : mCont (c) {}
    bool operator() (const typename Container::value_type &v) const
        { return (mCont.end() != std::find(mCont.begin(), mCont.end(), v)); }

private:
    Container mCont;
};

//------------------------------------------------------------------------------
/** Convenience function that allows construction of a recognize_t
    predicate without the need to specify template parameters. */
//------------------------------------------------------------------------------
template<class Container>
inline recognize_t<Container> recognize (const Container &c)
{
    return recognize_t<Container>(c);
}

//------------------------------------------------------------------------------
/** This allows string literals to be passed to recognize and is for
    backwards compatibility. */
//------------------------------------------------------------------------------
inline recognize_t<std::string> recognize (const char *c)
{
    return recognize_t<std::string>(c);
}

//------------------------------------------------------------------------------
/** Takes a sequence, an output iterator and a predicate, and breaks
    the sequence into smaller sequences according to the separators
    which are defined by the provided predicate.  The predicate should
    evaluate to true when applied to a separator value that is
    contained in the sequence. */
//------------------------------------------------------------------------------
template <class InIt, class OutIt, class Pred>
void tokenize (InIt first, InIt last, OutIt result, Pred is_sep)
{
    InIt i = first;
    InIt token_end = first;

    do
    {
        // Eat seperators
        while (i != last && is_sep(*i))
            ++i;

        // Find next token
        token_end = std::find_if(i, last, is_sep);

        // Append token to result
        if (i != token_end)
        {
            typedef typename OutIt::container_type::value_type value_type;
            *result++ = value_type(i, token_end);
        }

        i = token_end;
    } while (i != last);
}

//------------------------------------------------------------------------------
/** Tokenizes a sequence.  Same as above, but uses is_space() as
    default predicate. */
//------------------------------------------------------------------------------
template <class InIt, class OutIt>
inline void tokenize (InIt first, InIt last, OutIt result)
{
    tokenize(first, last, result, is_space());
}

//------------------------------------------------------------------------------
/** Tokenizes a string.  This is for backwards compatibility, but the
    above forms should be preferred, since they do not require an
    extra copy. */
//------------------------------------------------------------------------------
template <class Pred>
inline std::vector<std::string> tokenize (const std::string &s, Pred p)
{
    std::vector<std::string> result;
    tokenize(s.begin(), s.end(), std::back_inserter(result), p);

    return result;
}

//------------------------------------------------------------------------------
/** Tokenizes a string.  This is for backwards compatibility.  Same as
    above, but uses is_space() as default predicate.  For some reason,
    g++ won't recognize a default parameter. */
//------------------------------------------------------------------------------
inline std::vector<std::string> tokenize (const std::string &s)
{
    return tokenize(s, is_space());
}

#endif
