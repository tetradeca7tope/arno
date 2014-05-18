c============================================================
c     First Segment:
c
      program fred   ! Program name - fairly arbitrary
c
c     An example of a simple program.
c
      implicit none
      integer i, j, k
      integer addup
c
      i = 1
      j = 2
      k = addup (i, j)   ! Function call
      call output (k)    ! Subroutine call
      end
c============================================================
c     Second Segment:
c
      integer function addup (a, b)
c
c     A simple function to add together two numbers.
c     NOT the usual way of doing it!
c
      implicit none
      integer a, b
c
      addup = a + b
      return
      end
c============================================================
c     Third Segment:
c
      subroutine output (d)
c
c     Simple subroutine to print a number.
c     Quite unnecessary to do it this way!
c
      implicit none
      integer d
c
      print *, 'Here is the answer: ', d
      return
      end
c============================================================

