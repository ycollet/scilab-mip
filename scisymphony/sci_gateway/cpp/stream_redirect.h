/* -*- c++ -*- (enables emacs c++ mode) */
/*========================================================================

 Copyright (C) 2009-2011 Yann Collette

 This file is a part of GETFEM++

 Getfem++ is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as
 published by the Free Software Foundation; either version 2.1 of the
 License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.
 You should have received a copy of the GNU Lesser General Public
 License along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301,
 USA.

 As a special exception, you may use this file as part of a free software
 library without restriction.  Specifically, if other files instantiate
 templates or use macros or inline functions from this file, or you compile
 this file and link it with other files to produce an executable, this
 file does not by itself cause the resulting executable to be covered by
 the GNU General Public License.  This exception does not however
 invalidate any other reasons why the executable file might be covered by
 the GNU General Public License.

 ========================================================================*/

#ifndef STREAM_REDIRECT_H
#define STREAM_REDIRECT_H

#include <sciprint.h>

#include <iostream>
#include <streambuf>
#include <string>

//////////////////////////
// For cout redirection //
//////////////////////////

class ScilabStream : public std::basic_streambuf<char>
{
public:
  ScilabStream(std::ostream &stream) : m_stream(stream)
  {
    m_old_buf = stream.rdbuf();
    stream.rdbuf(this);
  }
  ~ScilabStream()
  {
    // output anything that is left
    if (!m_string.empty())
      sciprint("symphony: %s\n",m_string.c_str());

    m_stream.rdbuf(m_old_buf);
  }

protected:
  virtual int_type overflow(int_type v)
  {
    if (v == EOF)
      {
	sciprint("symphony: %s\n",m_string.c_str());
	m_string.clear();
      }
    else
      m_string.push_back(v);
    
    return v;
  }
  
  virtual std::streamsize xsputn(const char *p, std::streamsize n) 
  {
    m_string.append(p, p + n);
    
    int pos = 0;
    while (pos != std::string::npos)
      {
	pos = m_string.find(EOF);
	if (pos != std::string::npos)
	  {
	    std::string tmp(m_string.begin(), m_string.begin() + pos);
	    sciprint("symphony: %s\n",tmp.c_str());
	    m_string.erase(m_string.begin(), m_string.begin() + pos + 1);
	  }
      }
    
    return n;
  }
  
private:
  std::ostream   &m_stream;
  std::streambuf *m_old_buf;
  std::string     m_string;
};
#endif
