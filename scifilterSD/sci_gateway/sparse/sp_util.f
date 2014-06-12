christen this file util.f

c  Copyright (C) 1996 Roger Fletcher

c  Current version dated 06 January 1998

c  THE ACCOMPANYING PROGRAM IS PROVIDED UNDER THE TERMS OF THE ECLIPSE PUBLIC
c  LICENSE ("AGREEMENT"). ANY USE, REPRODUCTION OR DISTRIBUTION OF THE PROGRAM
c  CONSTITUTES RECIPIENT'S ACCEPTANCE OF THIS AGREEMENT

c*********************** dense matrix utilities ************************

      subroutine sp_rsol(n,nn,nmax,R,b)
      implicit double precision (a-h,o-z)
      dimension R(*),b(*)
c  solves Rx=b where R is nxn upper triangular. Solution overwrites b.
c  R is a single suffix array: the first nmax elements contain the first row
c  of R in positions 1:n, the next nmax-1 elements contain the second row of R,
c  and so on. nn indexes the element R(n,n) (where nn=n*(3-n)/2+(n-1)*nmax)
      n1=nmax+1
      ii=nn
      b(n)=b(n)/R(nn)
      do i=n-1,1,-1
        ii=ii-n1+i
        b(i)=-sp_scpr(-b(i),R(ii+1),b(i+1),n-i)/R(ii)
      enddo
      return
      end

      subroutine sp_rtsol(n,nn,nmax,R,b)
      implicit double precision (a-h,o-z)
      dimension R(*),b(*)
c  solves Rt.x=b with same conventions as above
c  nn is not required on entry but is set on exit
      n2=nmax+2
      nn=1
      b(1)=b(1)/R(1)
      do i=2,n
        i1=i-1
        call sp_mysaxpy(-b(i1),R(nn+1),b(i),n-i1)
        nn=nn+n2-i
        b(i)=b(i)/R(nn)
      enddo
      return
      end

      subroutine sp_Qprod(n,nmax,Q,x,b)
      implicit double precision (a-h,o-z)
      dimension Q(*),x(*),b(*)
c  forms b=M.x where Q is nxn, stored by columns, with stride nmax
      do i=1,n
        b(i)=0.D0
      enddo
      i1=1
      do i=1,n
        call sp_mysaxpy(x(i),Q(i1),b,n)
        i1=i1+nmax
      enddo
      return
      end

      subroutine sp_Qtprod(n,nmax,Q,x,b)
      implicit double precision (a-h,o-z)
      dimension Q(*),x(*),b(*)
c  forms b=M'.x where Q is nxn, stored by columns, with stride nmax
      i1=1
      do i=1,n
        b(i)=sp_scpr(0.D0,Q(i1),x,n)
        i1=i1+nmax
      enddo
      return
      end

      subroutine sp_brots(n,nmax,k,kk,R,v)
      implicit double precision (a-h,o-z)
      dimension R(*),v(*)
      ipip=kk
      do i=k-1,1,-1
        ip=i+1
        ipi=ipip-nmax+i
        ii=ipi-1
        call sp_angle(v(i),v(ip),cos,sin)
        call sp_rot(n-i,R(ipi),R(ipip),cos,sin)
        v(ip)=sin*R(ii)
        R(ii)=cos*R(ii)
        ipip=ii
      enddo
      return
      end

      subroutine sp_frots(nr,nc,nmax,R,v)
      implicit double precision (a-h,o-z)
      dimension R(*),v(*)
c nr is either nc or nc+1
      ii=1
      do i=1,nc
        ip=i+1
        ipi=ii+1
        ipip=ipi+nmax-i
        call sp_angle(R(ii),v(ip),cos,sin)
        call sp_rot(nr-i,R(ipi),R(ipip),cos,sin)
        ii=ipip
      enddo
      return
      end

      subroutine sp_angle(a,b,cos,sin)
      implicit double precision (a-h,o-z)
      z=sqrt(a**2+b**2)
      if(z.eq.0.D0)then
        cos=1.D0
        sin=0.D0
        return
      endif
      cos=a/z
      sin=b/z
      a=z
      b=0.D0
      return
      end

      subroutine sp_rot(n,a,b,cos,sin)
      implicit double precision (a-h,o-z)
      dimension a(*),b(*)
      if(sin.eq.0.D0)then
        if(cos.gt.0.D0)then
          do i=1,n
            b(i)=-b(i)
          enddo
        else
          do i=1,n
            a(i)=-a(i)
          enddo
        endif
      elseif(cos.eq.0.D0)then
        if(sin.ge.0.D0)then
          do i=1,n
            z=a(i)
            a(i)=b(i)
            b(i)=z
          enddo
        else
          do i=1,n
            z=a(i)
            a(i)=-b(i)
            b(i)=-z
          enddo
        endif
      else
        do i=1,n
          z=a(i)
          a(i)=cos*z+sin*b(i)
          b(i)=sin*z-cos*b(i)
        enddo
      endif
      return
      end

      subroutine sp_mysaxpy(a,x,y,n)
      implicit double precision (a-h,o-z)
      dimension x(*),y(*)
      if(a.eq.0.D0)return
      do i=1,n
        y(i)=y(i)+a*x(i)
      enddo
      return
      end

      subroutine sp_saxpys(a,x,is,y,n)
      implicit double precision (a-h,o-z)
c  saxpy with stride
      dimension x(*),y(*)
      if(a.eq.0.D0)return
      ix=1
      do i=1,n
        y(i)=y(i)+a*x(ix)
        ix=ix+is
      enddo
      return
      end

      subroutine sp_saxpyx(a,x,y,n)
      implicit double precision (a-h,o-z)
c  saxpy with result in x
      dimension x(*),y(*)
      if(a.eq.0.D0)then
        do i=1,n
          x(i)=y(i)
        enddo
      else
        do i=1,n
          x(i)=y(i)+a*x(i)
        enddo
      endif
      return
      end

      subroutine sp_saxpyz(a,x,y,z,n)
      implicit double precision (a-h,o-z)
c  saxpy with result in z
      dimension x(*),y(*),z(*)
      if(a.eq.0.D0)then
        do i=1,n
          z(i)=y(i)
        enddo
      else
        do i=1,n
          z(i)=y(i)+a*x(i)
        enddo
      endif
      return
      end

      subroutine sp_saxpyi(a,x,y,n)
      implicit double precision (a-h,o-z)
c  saxpy with interchange of x and y
      dimension x(*),y(*)
      if(a.eq.0.D0)then
        do i=1,n
          call sp_rexch(x(i),y(i))
        enddo
      else
        do i=1,n
          z=y(i)
          y(i)=x(i)+a*y(i)
          x(i)=z
        enddo
      endif
      return
      end

      function sp_scpr(a,x,y,n)
      implicit double precision (a-h,o-z)
      dimension x(*),y(*)
      sp_scpr=a
      do i=1,n
        sp_scpr=sp_scpr+x(i)*y(i)
      enddo
      return
      end

c     function sp_xlen(a,x,n)
c     implicit double precision (a-h,o-z)
c     dimension x(*)
c  finds the l_2 length of [a:x] where a is either 0.D0 or 1.D0
c  if overflow occurs the function is calculated in a less efficient way.
c  Users who cannot trap overflow should either use this method of calculation,
c  or use the alternative routine "xlen" below which is not quite so well
c  protected against overflow.
c     external  ieee_handler, abort
c     integer   ieee_flags, ieeer, ieee_handler
c     external  ieee_flags
c     character out*16
c     out = ''
c     ieeer = ieee_flags ( 'clearall','all','',out )
c     ieeer=ieee_handler('clear','overflow',abort)
c  this call of ieee_handler assumes that 
c         ieeer=ieee_handler('set','overflow',abort)
c  has been set in the driver. If not this call of ieee_handler and that below
c  should be removed
c     sp_xlen=a
c     do i=1,n
c       sp_xlen=sp_xlen+x(i)**2
c     enddo
c     sp_xlen=sqrt(sp_xlen)
c     ieeer=ieee_flags ( 'get','exception','',out )
c     if(out.eq.'overflow')then
c       call sp_linf(n,x,xmx,i)
c       xmx=max(xmx,1.D0) %this is needed if normalization is always used
c       sp_xlen=(a/xmx)**2
c       do i=1,n
c         sp_xlen=sp_xlen+(x(i)/xmx)**2
c       enddo
c       sp_xlen=xmx*sqrt(sp_xlen)
c       ieeer=ieee_flags ( 'clear','overflow','',out )
c     endif
c     ieeer=ieee_handler('set','overflow',abort)
c     return
c     end
      
      function sp_xlen(a,x,n)
      implicit double precision (a-h,o-z)
      dimension x(*)
      sp_xlen=a
      do i=1,n
        sp_xlen=sp_xlen+x(i)**2
      enddo
      sp_xlen=sqrt(sp_xlen)
      return
      end
      
      subroutine sp_linf(n,x,z,iz)
      implicit double precision (a-h,o-z)
      dimension x(*)
      z=0.D0
      do i=1,n
        a=abs(x(i))
        if(a.gt.z)then
          z=a
          iz=i
        endif
      enddo
      return
      end

      subroutine sp_r_shift(r,n,k)
      implicit double precision (a-h,o-z)
      dimension r(*)
      if(k.gt.0)then
        do i=1,n
          r(i)=r(i+k)
        enddo
      elseif(k.lt.0)then
        do i=n,1,-1
          r(i)=r(i+k)
        enddo
      endif
      return
      end

      subroutine sp_ishift(l,n,k)
      implicit double precision (a-h,o-z)
      dimension l(*)
      if(k.gt.0)then
        do i=1,n
          l(i)=l(i+k)
        enddo
      elseif(k.lt.0)then
        do i=n,1,-1
          l(i)=l(i+k)
        enddo
      endif
      return
      end

      subroutine sp_rexch(a,b)
      double precision a,b,z
      z=a
      a=b
      b=z
      return
      end

      subroutine sp_vexch(a,b,n)
      double precision a,b,z
      dimension a(*),b(*)
      do i=1,n
        z=a(i)
        a(i)=b(i)
        b(i)=z
      enddo
      return
      end

      subroutine sp_iexch(i,j)
      k=i
      i=j
      j=k
      return
      end
