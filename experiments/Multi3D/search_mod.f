c
c **********************************************************************
c
      logical function search(arr,x,ia,ib,k,err)
c
c * searches the point x within an error bound -err- in array arr      *
c * if that point is not found, it gives the                           *
c * interval in the array, where that point fits.                      *
c * the interval k if x in [arr(k),arr(k+1)), interval 1 if x in       *
c * (-inf,arr(2)), interval n-1 if x in [arr(n-1),inf)-(arr increasing)*
c * interval 1 if x in (inf,arr(2)),interval n-1 if x in [arr(n-1),    *
c * -inf)-(arr decreasing)                                             *
c * in the worst case a search costs 2*log2(n) cycles.                 *
c * this routine is very economic, when in a sequence of searches,     *
c * the differences between the points which are to be searched are    *
c * small.                                                             *
c * for example: a sequence of searches which goes along the array,    *
c * or a sequence of searches which remains in a small area of         *
c * the array.                                                         *
c *                                                                    *
c * input:     - arr         array with a non decreasing or a non      *
c *                          increasing fucntion.                      *
c *            - ia,ib       lower and upper limit of the array within *
c *                          is searched.                              *
c *            - x           value which you want to search            *
c *            - k           last index which was found                *
c *                          search                                    *
c *            - err         error bound                               *
c *                                                                    *
c * output:    - k           index of the point or interval of the     *
c *                          array which is found.                     *
c *            - search      true - x is found within the errorbound   *
c *                                  -err-                             *
c *                          false - x is not found.                   *
c *                                                                    *
c * send in by e.v.d.zalm, utrecht, sterrewacht, march 24th 1981       *
c **********************************************************************
      implicit real(8) (A-H,O-Z)
      dimension arr(ib)
c
      one=1.0
      search = .true.
      if (arr(ib)-arr(ia).lt.0.0) goto 100
5     if (x.lt.arr(ia).or.x.gt.arr(ib)) goto 25
      if (k.lt.ia.or.k.gt.ib-1) goto 25
      l = k
      ip = sign(one,x-arr(k))
      if (ip.lt.0.) goto 20
10    l = min(k+ip,ib)
      if (arr(l).ge.x) goto 30
      ip = ip * 2
      k = l            
      goto 10
20    k = max(l+ip,ia)
      if (arr(k).le.x) goto 30
      ip = ip * 2
      l = k
      goto 20
25    k = ia
      l = ib
30    if ((l-k).le.1) goto 50
      i = int((k+l)/2.)
      if (arr(i).ge.x) goto 40
      k = i
      goto 30
40    l = i
      goto 30
100   if (x.lt.arr(ib).or.x.gt.arr(ia)) goto 250
      if (k.lt.ia.or.k.gt.ib-1) goto 250
      l = k
      ip = sign(one,arr(k)-x)
      if (ip.lt.0.) goto 200
110   l = min(k+ip,ib)
      if (arr(l).le.x) goto 300
      ip = ip * 2
      k = l
      goto 110
200   k = max(l+ip,ia)
      if (arr(k).ge.x) goto 300
      ip = ip * 2
      l = k
      goto 200
250   k = ia
      l = ib
300   if ((l-k).le.1) goto 50
      i = int((k+l)/2.)
      if (arr(i).le.x) goto 400
      k = i
      goto 300
400   l = i
      goto 300
50    if (abs(x-arr(k)).le.err) return
      if (abs(x-arr(l)).gt.err) goto 60
      k = l
      return
60    search = .false.
      return
      end
c
c***********************************************************************
c