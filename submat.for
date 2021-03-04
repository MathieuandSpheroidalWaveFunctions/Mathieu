      module param
      integer, parameter :: knd = selected_real_kind(8)
      end module param
c      
         subroutine matfcn(lnum,ioprad,izxi,icq,isq,qc,r,iopang,narg,
     1                     arg,mc1c,mc1e,mc1dc,mc1de,mc23c,mc23e,mc23dc,
     2                     mc23de,ms1c,ms1e,ms1dc,ms1de,ms23c,ms23e,
     3                     ms23dc,ms23de,ce,ced,se,sed)

c      version 1.07 January 2021
c
c  Subroutine version of the fortran program matfcn originally developed
c  about 2005 by arnie lee van buren and jeffrey boisvert. Updated
c  several times since then. For more information see the GitHub
c  repository: GitHub.com/MathieuandSpheroidalWaveFunctions/Mathieu
c
c  purpose    : To calculate the Mathieu radial functions and their
c               first derivatives for a range of lnum orders from l
c               = 0 to l = lnum - 1 and for a given value of q
c               [or c=sqrt(2*q)] and the radial coordinate r. The
c               radial coordinate can either be the traditional
c               coordinate z or the spheroidal-like radial
c               coordinate xi = cosh(z).
c               To calculate the Mathieu angular functions and
c               their first derivatives for orders l = 0 to l = lnum - 1
c               for a given value of q (or c) and for one or more
c               angular coordinate values.
c
c  Matfcn is written in fortran 95. It is designed around the maximum
c  number of decimal digits ndec and the maximum exponent nex available
c  in real arithmetic. Procedures used in matfcn allow for exponents much
c  larger than nex since the resulting floating point radial function
c  values are given as a characteristic and an integer exponent. Testing
c  of matfcn shows that in real*8 arithmetic,it provides good results
c  for q positive for virtually all values of xi and for c up to at
c  least 10000. It also provides good results for q negative for
c  virtually all values of xi down to 1.000005 and for c up to at
c  least 10000i. See the discussion of accuracy below.
c
c  Matfcn can be run in either double precision or quadruple precision
c  arithmetic. The choice is set in the module param located at the
c  beginning of matfcn. Here, the kind parameter knd is set by the
c  statement:
c      integer, parameter :: knd = selected_real_kind(8)
c  Set the value of knd in the parenthesis to either 8 for 64 bit
c  arithmetic (double precision) or to 16 for 128 bit arithmetic
c  (quadruple precision). Using quadruple precision will provide
c  higher accuracy but will considerably increase the run time.
c
c  An alternative to using module param is to create two versions of
c  matfcn, one for double precision and one for quadruple precision. To
c  do this, replace knd throughout matfcn with the appropriate value of
c  kind for the desired precision. Also remove all occurrences of the
c  statement 'use param'.
c
c  Some computers may have more than 8 bytes for double precision
c  data and more than 16 bytes for quadruple precision data. In this
c  case just use the appropriate integers for the kind parameters in
c  module param. Also change the values of kindd and kindq set in
c  statement 5 below below the comments section to the number of
c  bytes for double precision data and quadruple precision data,
c  respectively. Larger values of kindd and kindq will lead to
c  more accurate results, especially for large values of c and small
c  values of x. Note that the mathematical constant gamma is given below
c  to 39 decimal digits. Of course, this will be truncated to the
c  correct precision for the kind parameter used. If matfcn is run
c  with higher precision than this, additional digits will be necessary
c  for gamma in order to ensure the higher accuracy associated with
c  this higher precision.
c
c  Some computer also use values for kind that do not correspond to the
c  number of bytes used in real data, In this case just use the values
c  for kind that correspond to double and quadruple precison arithmetic.
c  This includes setting kindd and kindq to the proper values for double
c  and quadruple precision arithmetic, respectively.
c
c    Input and output parameters appearing in the subroutine call
c    statement are defined below:
c
c          lnum   : number of integer values of l, given by
c                   0, 1, 2, ..., lnum-1, that radial and/or
c                   angular functions are desired
c
c          ioprad : (integer)
c                 : =0 if radial functions are not desired
c                 : =1 if radial functions of both kinds and
c                      their first derivatives are computed.
c                      When q is positive, this results in
c                      radial functions of the first kind Mc1 and
c                      Ms1 and of the second kind Mc2 and Ms2 and
c                      their first derivatives. When q is negative,
c                      the radial functions of the second kind are
c                      replaced by those of the third kind Mc3
c                      and Ms3 and their first derivatives.
c
c          izxi   : (integer)
c                 : = 1 if the input radial coordinate is the
c                       traditional radial coordinate z
c                 : = 2 if the input radial coordinate is the
c                       spheroidal-like radial coordinate = cosh(z)
c                       (actually x1 = xi - 1 is imput to avoid
c                       subtraction errors in the calculation of
c                       xi*xi - 1 in matfcn when xi is near unity.
c                       We avoid the subtraction error by using
c                       instead x1*(x1 + 2). Note that When z is
c                       input, the subtraction error is avoided in
c                       computing x1 = xi - 1 by using the identity
c                       x1=2*sinh(z/2)*sinh(z/2)
c
c          icq    : integer equal to 1 when q is input and equal to 2
c                   when cm, the magnitude of c, is input
c                   Note that c is pure imaginary when q is negative;
c                   we choose the positive square root of q in the
c                   definition of c.
c
c          isq    : integer equal to +1 when q is positive and
c                   equal to -1 when q is negative
c
c          qc     : size parameter [real(knd)]
c                 : = q when icq = 1
c                 : = cm, magnitude of c, when icq = 2
c                   This is equal to c when isq = +1 so that q
c                   is positive. When isq = -1 so that q
c                   is negative, then c is equal to i times the
c                   input value
c
c          r      : radial coordinate [real(knd)]
c                 : = z [real(knd)] when izxi = 1; this is the
c                   traditional radial coordinate ranging from
c                   zero to infinity.
c                 : = xi - 1 [real(knd)] when izxi = 2; where xi is
c                   the spheroidal-like radial coordinate equal to
c                   cosh(z). A dummy value, e.g. 1.0d0, should be
c                   entered for r if ioprad = 0
c
c          iopang : (integer)
c                 : = 0 if angular functions are not desired
c                 : = 1 if angular functions of the first kind
c                       are computed
c                 : = 2 if angular functions of the first kind and
c                       their first derivatives are computed
c
c          narg   : number of values of the angular coordinate phi for
c                   which angular functions of the first kind are
c                   calculated (integer). If no angular functions
c                   are desired (iopang = 0), set narg equal to unity.
c
c          arg:     vector arg(narg) containing the values of phi in
c                   degrees for which angular functions are desired
c                   [real(knd)]
c
c          mc1c   : real(knd) vectors of length lnum containing the
c          mc1dc    characteristics for the cosine radial functions
c                   of the first kind Mc1 and their first derivatives
c                   with respect to z.
c                   When q is positive, the functions are real.
c                   When q is negative and the order is odd, the
c                   functions are imaginary. The real values given
c                   in the vectors dmc1 and dmc1d must be multiplied
c                   by i to obtain the characteristics for odd orders.
c
c          mc1e   : integer vectors of length lnum containing the
c          mc1de    exponents corresponding to mc1c and mc1dc.
c
c          mc23c  : When q is positive, these are real(knd) vectors of
c          mc23dc   length lnum containing the characteristics for the
c                   cosine radial functions of the second kind Mc2 and
c                   their first derivatives with respect to z.
c                   When q is negative, these vectors instead contain
c                   the characteristics for the cosine radial functions
c                   of the third kind Mc3 and their first derivatives
c                   with respect to z. Mc3 and its first derivative is
c                   real when l is odd and imaginary when l is even.
c                   In this case, the real values given in the vectors
c                   dmc23 and dmc23d must be multiplied by i to obtain
c                   the characteristics for even orders.
c
c          mc23e :  integer vectors of length lnum containing the
c          mc23de   exponents corresponding to dmc23 and dmc23d
c
c          ms1c   : real(knd) vectors of length lnum containing the
c          ms1dc    characteristics for the sine radial functions
c                   of the first kind Ms1 and their first derivatives
c                   with respect to z.
c                   When q is positive, the functions are real.
c                   When q is negative and the order is odd, the
c                   functions are imaginary. The real values given in
c                   the vectors dms1 and dms1d must be multiplied by i
c                   to obtain the characteristics for odd orders.
c
c          ms1e   : integer vectors of length lnum containing the
c          ms1de    exponents corresponding to ms1c and ms1de
c
c          ms23c  : When q is positive, these are real*8 vectors of
c          ms23dc   length lnum containing the characteristics for the
c                   sine radial functions of the second kind Ms2 and
c                   their first derivatives with respect to z.
c                   When q is negative, these vectors instead contain
c                   the characteristics for the sine radial functions
c                   of the third kind Ms3 and their first derivatives
c                   with respect to z. Ms3 and its first derivative is
c                   real when l is odd and imaginary when l is even.
c                   In this case, the real values given in the vectors
c                   dms23 and dms23d must be multiplied by i to obtain
c                   the characteristics for even orders.
c
c          mc23e  : integer vectors of length lnum containing the
c          mc23de   exponents corresponding to dms23 and dms23d
c
c          ce,ced : vectors ce(lnum,narg) and ced(lnum,narg) that
c                   contain narg calculated angular cosine functions
c                   and their first derivatives for each of the lnum
c                   values of l [real(knd)]
c                   For example, ce(10,1) is the angular function for
c                   l = 9 and the first value of the angle phi given
c                   by arg(1)
c          se,sed : vectors se(lnum,narg) and sed(lnum,narg) that
c                   contain narg calculated angular sine functions
c                   and their first derivatives for each of the lnum
c                   values of l [real(knd))
c                   For example, se(10,1) is the angular function for
c                   l = 9 and the first value of the angle phi given
c                   by arg(1)
c
c  We use the term radial function for what is normally referred to as
c  a modified Mathieu function. We do this to be consistent with our
c  usage of the term radial function for spheroidal wave functions.
c  The use of what we call the size parameter c instead of the
c  traditional q is consistent with usage for spheroidal functions
c  (see, e.g., Handbook of Mathematical Functions). So is the use of
c  xi = cosh(z) as the radial coordinate. We note that the
c  algebraic form of Mathieu's equation can be obtained from the
c  differential equation for radial prolate spheroidal functions by
c  setting m = 1/2.
c
c  The first derivative values provided for the radial functions
c  are with respect to the traditional radial coordinate z that ranges
c  from 0 to infinity. First derivatives with respect to the radial
c  coordinate xi can be obtained by multiplying those given by matfcn
c  by the factor 1/sqrt(xi*xi-1). See the comments below about using
c  x1 = xi - 1 as input to matfcn to avoid subtraction errors in forming
c  the expression xi*xi - 1 when xi is near unity. Here xi*xi - 1 is
c  then accurately computed using x1*(x1 + 2).
c
c  The first derivative values provided for the angular functions
c  are with respect to the traditional angle coordinate.
c
c  The methods used in matfcn are documented in the journal article:
c  A. L. Van Buren and J. E. Boisvert, "Accurate calculation of the
c  modified mathieu functions of integer order," Quart. Appl.
c  Math. 65, 1-23 (2007). A pdf version of this article is available
c  for download on the website indicated above. [Note that the argument
c  for both cos and sin in the two equations in (10) is in error and
c  should be r times theta.]
c
c          Accuracy of Results Using Real*8 Arithmetic
c
c  When c is real (q positive) and xi is not zero, matfcn provides
c  values for the radial functions of both the first and second kinds
c  that are usually accurate to at least 10 decimal digits. Matfcn
c  estimates the accuracy as the nummber of decimal digits of agreement
c  between the theoretical Wronskian and its value computed using the
c  resulting radial function values. We tested matfcn for a selection
c  of xi values ranging from 1.00000000001 to 100 and for a range of c
c  values up to 1000. We also tested all of the xi values up to 10 for
c  c values up to 10000. Lnum was chosen to be the larger of 1000 or an
c  order where the radial functions of the first kind had magnitudes
c  well below 10 to the -300 power. The resulting accuracies were 10 or
c  more digits with c up to 100. Accuracies were at least 8 digits for c
c  = 200 with only three occurrences of 8 digits. For c = 500, there
c  were three 7 digit results. For larger c there were a few 7 digit
c  and two 6 digit results. Matfcn is expected to provide useful results
c  for xi outside the test range and for lnum much higher than tested.
c  It is always possible to choose xi so close to a root that lower
c  accuracy is obtained. The corresponding radial function would then
c  have a somewhat smaller magnitude than neighboring function values
c  and likely would not adversely affect solutions to problems involving
c  these functions.
c
c  When c is real and xi is equal to unity, the sine radial functions
c  of the first kind and the first derivative of the cosine radial
c  function of the first kind are equal to zero. The other radial
c  functions of the first kind are very accurate unless extremely close
c  to a root. When c is large and the order l is somewhat less than
c  2c/pi, both the cosine radial functions of the second kind and the
c  first derivatives of the sine radial functions of the second kind
c  tend to have a very small magnitude. For a given large value of c,
c  the magnitude decreases as the order decreases toward zero. For a
c  given order, the magnitude decreases as c increases. These functions
c  are inaccurate to the degree that they are small. Thus they should
c  not reduce the accuracy of solutions to problems where the radial
c  functions of the second kind appear in a pair with the functions of
c  the first kind. Here the accurate radial functions of the first kind
c  are larger and dominate. Note that the Wronskian accuracy check
c  overestimates the accuracy here.
c
c  When c is imaginary (q negative) and xi is not equal to zero, matfcn
c  provides values for the radial functions of the first and third kinds
c  that are usually accurate to 10 or more digits. We tested matfcn
c  using a set of xi values ranging from 1.000000000001 to 100 and a set
c  of c values up to 10000i. We tested the full set of xi values for c
c  up to 10i. We tested the xi values ranging from 1.000005 to 100 for c
c  up to 5000i and the xi values ranging from 1.000005 to 10 for c up to
c  10000i. Lnum was chosen to be the larger of 1000 or a value large
c  enough that the radial functions of the first kind have magnitudes
c  well below 10 to the -300 power. The resulting accuracies were 10 or
c  more digits for c up to 50i and 9 or more digits for c up to 500i.
c  For higher values of c, accuracies were 8 or more digits except for a
c  few 7 digit results and one 6 digit result.
c
c  For c larger than about 10i, matfcn uses the Bessel function (K)
c  expansion to compute radial functions of the third kind for orders
c  less than 2*cm/pi. The Bessel product expansion suffers too much
c  subtraction error here to provide accurate results. The number of
c  terms maxk required for convergence of the K expansion to the
c  desired accuracy and the time to perform the calculations increases
c  somewhat linearly with 1/(xi-1). For xi = 1.000001 and c = 10000i,
c  storing the vectors of expansion functions and coefficients for the
c  calculations can require as much as 1 GB of memory. I set the maximum
c  value of maxk to 10000000 (maxterm) for my testing. The user can
c  increase this in statement number 20 below if they desire smaller xi
c  values for c greater than 10i and have sufficient computer memory.
c  Note that the necessary memory doubles for real*16 arithmetic and
c  that the Bessel function expansion can be used to higher values
c  of c before the K expansion is required.
c
c  When c is imaginary and xi is equal to unity, the cosine radial
c  function of the first kind and the first derivative of the sine
c  radial function are equal to zero. The Bessel function expansion
c  used for radial functions of the third kind for orders below about
c  2(cm)/pi is useless for xi = 0. Unfortunately the Bessel product
c  expansion that must be used here suffers subtraction errors that
c  are too large to allow accurate values for low order radial functions
c  of the third kind when c is somewhat greater than 20i. Note
c  that when the function accuracy is such that the Wronskian accuracy
c  is less than the desired minimum accuracy, matfcn uses the
c  Wronskian to compute an accurate value for the sine radial
c  function of the third kind and the first derivative of the cosine
c  radial function of the third kind. Of course, the cosine radial
c  function and the first derivative of the sine radial function
c  likely have an accuracy less than the desired minimum accuracy and
c  in fact are set equal to zero when the Wronskian aaccuracy is zero.
c  Radial function values for orders near and above 2(cm)/pi are
c  usually highly accuracy for all values of c.
c
c  The calculated angular functions are highly accurate except for
c  orders less than about 2(cm)/pi when c is large in magnitude.
c  Function values tend to be less accurate the larger the magnitude of
c  c, the smaller the order and the closer the angle is to zero degrees.
c  However, the loss in accuracy is due to subtraction error in the
c  calculations and is accompanied by a proportional decrease in the
c  magnitude of the angular function value relative to the corresponding
c  trigonometric function values. This decrease in magnitude almost
c  always results in a corresponding reduction in the magnitude of their
c  contribution to the solution of physical problems involving Mathieu
c  functions. Note that the angular functions have the same norm pi as
c  the corresponding cosine and sine functions.
c                                                         
c          Accuracy of Results Using Real*16 Arithmetic
c
c  If the user desires more accuracy than is provided using real*8
c  arithmetic, matfcn can be run using real*16 arithmetic. However,
c  the execution time increases considerably. Matfcn computes results
c  to full precision except for the radial functions of the third kind
c  m3 for orders less than 2*cm/pi. These are computed to a precision
c  based on the parameter minacc. Minacc is set below after these
c  comments to be equal to 15 digits for real*16 arithmetic. It can be
c  increased to provide greater precision and accuracy for lower order
c  values of m3. Although when c is imaginary and large in magnitude and
c  xi is very near unity, the accuracy may not increase much when using
c  the expansion in K functions. Here, many more terms than used in
c  matfcn would be required to obtain significantly higher accuracy.
c  It is advised that one not change minacc from its value of 10 given
c  below when using real*8 arithmetic.
c
c  If desired, output data can be written to one or more of the
c  following five files: fort.20, fort.30, fort.40, fort.50 and fort.60.
c  The write and corresponding format statements are commented out
c  below to allow users to avoid these files. The user can obtain any
c  desired output file, say fort.20, by searching matfcn for write(20,
c  and the associated format statements and removing the c in column 1
c  of the lines of code found. Some compilers also require that the file
c  be opened. In that case remove the c in the line of code below
c  opening the file. Users running matfcn outside the range over which
c  it was tested, especially when the size parameter c has a magnitude
c  larger than 10000 may want to consider using fort.60. Matfcn writes
c  to fort.60 the values for xi, c, m, l, and estimated accuracy when
c  the estimated accuracy falls below a specified number of integer
c  digits. See comments below about this file.
c
c   fort.20
c
c     This file contains values for all radial functions that have
c     been calculated. The first line contains the values for either
c     z or xi, depending on whether z or xi - 1 is input and either c
c     (the factor i = sqrt(-1) is suppressed from c when q is negative)
c     or q, depending on whether c or q is input, formatted as follows
c     (see statements 20 and 30 in the first subroutine main):
c
c                z or xi : e24.15 for real*8; e40.31 for real*16
c                c or q  : e24.15 for real*8; e40.31 for real*16
c
c     Each subsequent line in matfcn.rad contains radial functions
c     for given values of l, beginning with 0 and continuing to l =
c     lnum - 1. First is given the value of l, followed on the
c     same line first by the cosine radial functions and then by a
c     two digit estimate of the number of accurate digits in the
c     preceeding radial functions. For l unequal to zero, the line
c     of cosine radial functions is followed by a line of the
c     corresponding sine radial functions for that value of l
c     together with a two digit estimate of their accuracy.
c     the estimate of accuracy is determined by comparing the
c     theoretical Wronskian to the value computed using the radial
c     functions and their first derivatives. The output and
c     corresponding format for each line is as follows
c     (see statements 520, 530, 550, and 560 in subroutine main)
c
c        for ioprad = 1
c
c               l      : value for l (i5)
c               mc1c/  : characteristic of the cosine/sine radial
c               ms1c     function of first kind (f17.14). when q is
c                        negative (c is imaginary) and l is odd, mc1c
c                        and ms1c are pure imaginary and the factor i
c                        is suppressed
c               imc1e/ : exponent of the cosine/sine radial function of
c               ims1e    first kind (i6 for q positive; i7 for q
c                        negative)
c               mc1dc/ : characteristic of the first derivative of the
c               ms1dc    radial function of first kind (f17.14). when q
c                        is negative and l is odd, mc1dc and ms1dc are
c                        pure imaginary and the factor i is suppressed
c               imc1e/ : exponent of the first derivative of the
c               ims1e    radial function of first kind (i6 for q
c                        positive; i8 for q negative)
c               mc23c/ : characteristic of the radial function
c               ms23c    of second (or third) kind (f17.14). when q
c                        is negative and l is even, mc3c and ms3c are
c                        pure imaginary and the factor i is suppressed
c               imc23e/: exponent of the radial function of the
c               ims23e   second (or third) kind (i6 for q positive;
c                        i8 for q negative)
c               mc23dc/: characteristic of the first derivative of the
c               ms23dc   radial function of the second (or third) kind
c                        (f17.14). when q is negative and l is even,
c                        mc3dc and ms3dc are pure imaginary and the
c                        factor i is suppressed
c               imc23de/:exponent of the first derivative of the
c               ims23de  radial function of the second (or third)
c                        kind (i6 for q positive; i8 for q negative)
c               naccr  : accuracy: equal to the number of decimal digits
c                        of agreement between the theoretical Wronskian
c                        and the calculated Wronskian (i2)
c
c     When one of the two products of radial functions used in forming
c     the Wronskian is significantly smaller in magnitude than the other
c     product, the Wronskian accuracy naccr might overestimate the
c     accuracy of one or both of the two functions in the smaller
c     product. They can be less accurate than naccr by as many digits
c     as obtained by integer truncation of the logarithm to the base 10
c     of the ratio of the magnitude of the larger product to the smaller
c     product. The file fort.40 (described below) will show the
c     subtraction errors involved in the calculation of the functions
c     and identify inaccuracies in their values.
c
c   fort.30
c
c     This file contains values for all angular functions that have
c     been calculated. The first line contains the value for either c
c     (i is suppressed when q is negative) or q, depending on which
c     one is input, formatted as follows (see statements 40 and 45
c     in subroutine main):
c
c                c or q : e24.15 in real*8; e40.31 for real*16
c
c     The second line contains the value for the first l (=0),
c     formatted as follows (see statement 90 in subroutine main):
c
c                l      : i5
c
c     This is followed by a series of narg lines. each line contains
c     a desired value of angle followed by the corresponding angular
c     functions and accuracy. Specific output and format for each line
c     is as follows:
c
c        for iopang = 1:
c
c               arg    : for ioparg = 0, angle in degrees (f20.14; see
c                        statements 620 and 630 in subroutine main)
c               ce     : value of the angular function ce (e24.15;
c                        see statements 620 and 630 in subroutine main)
c               se     : value of the angular function se (e24.15;
c                        see statements 620 and 630 in subroutine main).
c                        for l = 0 the value for se is set equal to zero.
c
c        for iopang = 2, each line also includes:
c               ced    : value of the first derivative of the
c                        angular function ce, listed immediately
c                        following the value for ce: (e24.15;
c                        see statement 630 in subroutine main)
c               sed    : value of the first derivative of the
c                        angular function se, listed immediately
c                        following the value for se: (e24.15);
c                        see statement 630 in subroutine main).
c                        for l = 0 the value for sed is set equal to
c                        zero.
c
c        for iopang = 1 or 2:
c      nacca           : accuracy: smaller of naccc and naccs, where
c                        naccc is the number of decimal digits of
c                        accuracy estimated for the angular function
c                        ce and its first derivative when iopang = 2
c                        and where naccs is the corresponding estimate
c                        for se. naccc and naccs are conservative
c                        estimates given by ndec -2 - the larger of
c                        the subtraction error encountered in the
c                        series calculation of the corresponding
c                        angular function and its first derivative.
c                        When nacca is equal to 0, the corresponding
c                        angular functions and their first derivatives
c                        are set equal to zero. naccc and naccs are
c                        written with format i2; see statements 620
c                        and 630 in subroutine main)
c
c   fort.40 and fort.50
c
c     These files are diagnostic files. They contain information about
c     the specific techniques used and the numbers of terms required
c     for the radial function and angular function calculations,
c     respecively. They are annotated and should be self-explanatory.
c     If fort.40 is not desired, then search for (40, in the program
c     and comment out the lines of code found by placing a c in column
c     one. Similary if fort.50 or fort.60 is not desired search for
c     (50, or (60 and comment out the lines of code found. The
c     associated open and format statements can also be commented out.
c
c   fort.60
c
c     This file may be of interest to the user. It is recommended that
c     the user utilize this file, especially when using matfcn for very
c     large imaginary values of c with x near unity. Whenever the
c     estimated accuracy for the radial functions falls below a
c     designated integer value, the associated values of x, c (and the
c     corresponding value for q), and l are written to fort.60. This
c     integer is currently set equal to 6 in the write statement for
c     fort.60 found just before the line numbered 600 below in
c     subroutine main. The reason for choosing 6 is that it is expected
c     that 6 accurate decimal digits are sufficient for most
c     applications. The integer can be changed to any desired value in
c     the write(60,*) statement below.
c
c     In the unlikely case that the eigenvalue routine fails to properly
c     converge, the value of l for which this occurs will be written
c     to fort.60. Note that this has never been observed with matfcn.
c
c          Expansion coefficients A and B
c
c  The user may be interested in the expansion coefficients A and B for
c  each order l. Both are calculated as vectors containing ratios of
c  successive coefficients. The vector enra(j) contains the ratios
c  A(sub 2j)/A(sub 2j-2) for l even and A(sub 2j+1)/A(sub 2j-1) for l
c  odd. Enra is computed in convera and returned to subroutine main.
c  Similarly, the vector enrb(j) contains the ratios B(sub 2j)/
c  B(sub 2j-2) for even l and B(sub 2j+1)/B(sub 2j-1) for l odd. Enrb is
c  computed in converb and also returned to subroutine main. Note that
c  B(sub 0) is not used so the ratio enrb(1) for l even is also not
c  used.
c
c  The subroutines convera and converb return coefficients appropriate
c  for q positive. If q is negative, the A and B coeffficients are the
c  same as those for q positive when l is even. However, when l is odd,
c  the A coefficients for q negative are equal to the B coefficients for
c  q positive. And the B coefficients for q negative are equal to the A
c  coefficients for q positive. Thus the coefficient ratios enra and
c  enrb returned by convera and converb must be interchanged when q is
c  negative.
c
c  Both enra and enrb have the dimension maxd, which is the
c  estimate of the largest number of coefficients that will be required
c  to compute the radial and/or angular functions for all desired lnum
c  orders. The number of coeffficient ratios computed for each l is 
c  given by limd/2, where limd tends to increase with increasing l. It
c  can be very large when xi is very close to unity and the traditional
c  expansion in modified Bessel function of the second kind K is used to
c  compute the radial functions m3. If enra and enrb are to be used in
c  calculations involving just the the angular functions, then only
c  limdang/2 coefficients are needed for each l. The value of limdang is
c  determined for each l below between statements 90 and 95 in
c  subroutine main. Note that limdang is only determined when matfcn is
c  run with iopang equal to 1 or 2.
c
c  Matfcn computes the coefficient A(sub l) for order l in subroutine
c  cese where the angular functions are computed and returns it to
c  subroutine main as asubl. The A coefficient with the next higher
c  index A(sub l+2) is obtained by multiplying A(sub l) by
c  enra(sub (l+2)/2) for l even and by enra(sub (l+1)/2) for l odd.
c  Higher index cofficients are then obtained in order by multiplying by
c  successive ratios enra. The A coefficient A(sub l-2) is obtained by
c  dividing A(sub l) by the ratio enra(sub l/2) for even l and by
c  enra(sub (l-1)/2) for l odd. Lower order coefficients are then
c  obtained in turn by successive division by lower index ratios enra.
c  The B coefficients are obtained in a similar fashion starting with
c  B(sub l) computed in subroutine cese and returned to main as bsubl.
c  Note that the values for asubl and bsubl are correct as given for
c  both q positive and q negative. No interchange is needed when q is
c  negative and l is odd.  
c
c          Eigenvalues
c
c  The eigenvalues for the expansion coefficients A are computed in
c  subroutine convera and returned to main where they are stored in the
c  vector eiga(l+1). Similarly the eigenvalues for B are computed in
c  subroutine converb and returned to main and stored in the vector
c  eigb(l+1). Note that there is no B eigenvalue for l = 0 so eigb(1)
c  is set to zero here. Note also that for q negative the eigenvalues 
c  for A and B switch places when l is odd; thus the A eigenvalue is
c  computed in converb and the B eigenvalue is computed in convera. The
c  vectors eiga and eigb for q negative given below correctly reflect
c  this interchange so no further change needs to be made.
c
c  real(knd) scalars
c
        use param
c
        real(knd) qc,q,cm,r,z,x1
        real(knd) mc1c(lnum),mc1dc(lnum),mc23c(lnum),mc23dc(lnum),
     1            ms1c(lnum),ms1dc(lnum),ms23c(lnum),ms23dc(lnum),
     2            ce(lnum,narg),ced(lnum,narg),se(lnum,narg),
     3            sed(lnum,narg),arg(narg)
        integer   mc1e(lnum),mc1de(lnum),mc23e(lnum),mc23de(lnum),
     1            ms1e(lnum),ms1de(lnum),ms23e(lnum),ms23de(lnum)
c
c
c  ndec: the maximum number of decimal digits available in real(knd)
c           arithmetic.
c  nex:  the maximum exponent available in real(knd) arithmetic.
c
        ndec=precision(cm)
        nex=range(cm)-1
c
c  Here is where the user sets kindd, the number of bytes available
c  in double precision data for the computer that coblfcn is run on.
c  Similarly kindq, the number of bytes available in quadruple
c  precision data, is set here. This allows coblfcn to be run on
c  computers with values for kindd and kindq different than 8 and 16,
c  respectively.
c
5       kindd=8
        kindq=16

c  set the minimum desired accuray minacc to 10 for real*8
c  arithmetic and to 15 for real*16 arithmetic. Minacc for
c  real*16 can be changed if desired. See comments above about
c  changing minacc
c
        if(knd.eq.kindd) minacc=10
        if(knd.eq.kindq) minacc=15
c
c  open output files
c        open(20,file='fort.20')
c        open(30,file='fort.30')
c        open(40,file='fort.40')
c        open(50,file='fort.50')
        open(60,file='fort.60',access='append')
c
          if(icq.eq.1) then
          q=qc
          cm=sqrt(2.0e0_knd*abs(q))
          end if
          if(icq.eq.2) then
          cm=qc
          q=cm*cm/2.0e0_knd
          if(isq.eq.-1) q=-q
          end if
          if(izxi.eq.1) then
          z=r
          x1=2.0e0_knd*(sinh(0.5e0_knd*z))**2
          end if
          if(izxi.eq.2) then
          x1=r
          z=0.0e0_knd
          end if 
c
c  set array dimensions
        maxj=1
        if(ioprad.ne.0) maxj=1.5*lnum+4*ndec+int(cm)+105
c
c      for isq = -1, the maximum s value used for m1 occurs at the
c      last l that the Bessel product series is used, namely at the
c      break point l = lb = 2*cm/pi.
c
        if(isq.eq.1) go to 10
        ismax=int(0.85e0_knd*cm)
        if(x1.lt.1.0e0_knd) ismax=ismax*1.5
        if(x1.lt.0.1e0_knd.and.x1.ge.0.00001e0_knd) ismax=ismax*
     1                           int(3.0e0_knd**(-1.0e0_knd-log10(x1)))
        if(x1.lt.0.00001e0_knd) ismax=ismax*81
        maxj=max(maxj,ismax+int(cm/3.14159)+4*ndec+
     1           int(cm)+105)
10      maxp=2*lnum+2*cm+4*ndec+105
        maxn=maxj
        maxk=1
        maxkbp=1
        if(ioprad.eq.0.or.isq.eq.1.or.x1.eq.0.0e0_knd) go to 30
        maxkbp=lnum+4*ndec+int(cm)+103
        if(minacc.le.(ndec-2).and.cm.le.5.0e0_knd) go to 30
        if(minacc.le.(ndec-4).and.cm.le.10.0e0_knd) go to 30
        if(minacc.le.(ndec-6).and.cm.le.15.0e0_knd) go to 30
        if(minacc.le.(ndec-8).and.cm.le.20.0e0_knd) go to 30
        maxk=1.5*cm+100
        if(x1.lt.1.0e0_knd.and.x1.ge.0.01e0_knd) maxk=40/x1+cm/sqrt(x1)+
     1                                        1.5*cm+200
        if(x1.lt.0.01e0_knd.and.x1.gt.0.0e0_knd) maxk=34/x1+1.5*cm+
     1                                           1.4*cm/sqrt(x1)+200
20      maxterm=10000000
        maxk=min(maxk,maxterm)
        maxk=max(maxk,maxj)
30      maxd=max(maxj,maxn,maxk,maxkbp,maxp)/2+1
        maxlp=lnum+3
        if(isq.eq.-1) maxlp=max(maxlp,ismax+int(cm/3.14159)+3)
        ngau=200
c
        call main (lnum,cm,q,icq,isq,ioprad,iopang,minacc,izxi,x1,z,
     1             narg,arg,maxd,maxj,maxlp,maxn,maxp,maxkbp,maxk,
     2             ndec,nex,ngau,kindd,kindq,mc1c,mc1e,mc1dc,mc1de,
     3             mc23c,mc23e,mc23dc,mc23de,ms1c,ms1e,ms1dc,ms1de,
     4             ms23c,ms23e,ms23dc,ms23de,ce,ced,se,sed)
c
        end
c
c
        subroutine main (lnum,cm,q,icq,isq,ioprad,iopang,minacc,izxi,x1,
     1                   z,narg,arg,maxd,maxj,maxlp,maxn,maxp,maxkbp,
     2                   maxk,ndec,nex,ngau,kindd,kindq,amc1c,mc1e,
     3                   amc1dc,mc1de,amc23c,mc23e,amc23dc,mc23de,ams1c,
     4                   ms1e,ams1dc,ms1de,ams23c,ms23e,ams23dc,ms23de,
     5                   ace,aced,ase,ased)
c
c  purpose:    to coordinate the calculation of Mathieu functions,
c              including the radial and angular functions and their
c              first derivatives using various algorithms.
c
c     input:    lnum   : desired number of values of l = m, m + 1, ...,
c                        m + lnum - 1
c               cm     : magnitude of c, the spheroidal-like size
c                        parameter = sqrt(2*q)
c               q      : traditional size parameter
c               icq    : equal to 1 when q is input and equal to 2 when
c                        cm, the magnitude of c, is input
c               isq    : sign of q
c               ioprad : = 0 if radial functions are not desired
c                      : = 1 if radial functions of both kinds and
c                       their first derivatives are computed.
c                       When q is positive, this results in
c                       radial functions of the first and second
c                       kind m1 and m2. When q is negative, it
c                       results in radial functions of the first
c                       and third kind m1 and m3.
c               iopang : equal to 0 if no angular functions are desired;
c                        equal to 1 if only angular functions of the
c                          first kind are desired;
c                        equal to 2 if angular functions of the first
c                          kind and their first derivatives are desired
c               minacc : desired minimum accuracy for the radial
c                        functions
c               izxi   : = 1 if the input radial coordinate is the
c                       traditional radial coordinate z
c                      : = 2 if the input radial coordinate is the
c                       spheroidal-like radial coordinate xi = cosh(z)
c               x1     : radial coordinate xi - 1 = cosh(z) - 1
c               z      : traditional radial coordinate 
c               arg    : vector of narg angle coordinates in degrees
c                        for which angular functions are desired
c               narg   : number of desired angle arguments.                       
c               maxd   : dimension of the vectors enra and enrb that
c                        contain ratios of the a and b expansion
c                        coefficients
c               maxj   : equal to the dimension of the array of ratios
c                        of Bessel functions of the first kind and of
c                        the array of ratios of the first derivatives of
c                        these functions. The ratios are of consecutive
c                        functions of the same parity. Maxj applies to
c                        both the regular and the modified Bessel
c                        functions.
c               maxlp  : lnum+3
c               maxn   : dimension of the array of ratios of Bessel
c                        functions of the second kind and of the array
c                        of ratios of the first derivatives of these
c                        functions. The ratios are of consecutive
c                        functions of the same parity.                 
c               maxp   : number of sine and cosine function values
c                        computed for each desired angle coordinate
c               maxk   : dimension of the array of ratios of modified
c                        Bessel functions of the third kind and of the
c                        array of ratios of the first derivatives of
c                        these functions. The ratios are of consecutive
c                        functions of the same parity. 
c               maxkbp : dimension of the array of ratios of modified
c                        Bessel functions of the third kind and of the
c                        array of ratios of the first derivatives of
c                        these functions. The ratios are of consecutive
c                        functions of the same parity.                                               
c               ndec   : number of decimal digits for real(knd)
c               nex    : maximum exponent for real(knd)
c               ngau   : order of the Gaussian quadrature to be used in
c                        computing integrals in subroutine pint for use
c                        in subroutine r2int where the integal method
c                        is used to calculate r2 and r2d
c               kindd  : kind value for double precision real data
c               kindq  : kind value for quadruple precision real datac
c
c     Output:   amc1c  : vector containing the lnum characteristics for
c                        the cosine radial functions of the first kind
c                        Mc1. The factor i is suppressed for all of the
c                        radial and angular functions that are maginary. 
c               mc1e   : integer vector of exponents corresponding to
c                        amc1c
c               amc1dc : vector containing the lnum characteristics for
c                        the first derivatives with respect to z of the
c                        cosine radial functions of the first kind
c               mc1de  : integer vector of exponents corresponding to
c                        amc1dc                       
c               amc23c : When q is positive, these are vectors of lnum
c                        characteristics for the cosine radial functions
c                        of the second kind Mc2.
c                        When q is negative, these are vectors of lnum
c                        characteristics for the cosine radial functions
c                        of the third kind Mc3
c               mc23e  : integer vector of exponents corresponding to
c                        amc23c
c               amc23dc: When q is positive, this is a vector of lnum
c                        characteristics for the first derivative with
c                        respect to z of the cosine radial functions
c                        of the second kind
c                        When q is negative, this is a vector of lnum
c                        characteristics for the first derivatives with
c                        respect to z of the cosine radial functions
c                        of the third kind
c               mc23de : integer vector of exponents corresponding to
c                        amc23dc
c               ams1c  : vector containing the lnum characteristics for
c                        the sine radial functions of the first kind
c                        Ms1
c               ms1e   : integer vector of exponents corresponding to
c                        ams1c
c               ams1dc : vector containing the lnum characteristics for
c                        the first derivatives with respect to z of the
c                        sine radial functions of the first kind
c               ms1de  : integer vector of exponents corresponding to
c                        ams1dc                       
c               ams23c : When q is positive, these are vectors of lnum
c                        characteristics for the sine radial functions
c                        of the second kind Ms2.
c                        When q is negative, these are vectors of lnum
c                        characteristics for the sine radial functions
c                        of the third kind Ms3
c               ms23e  : integer vector of exponents corresponding to
c                        ams23c
c               ams23dc: When q is positive, this is a vector of lnum
c                        characteristics for the first derivative with
c                        respect to z of the sine radial functions
c                        of the second kind
c                        When q is negative, this is a vector of lnum
c                        characteristics for the first derivatives with
c                        respect to z of the sine radial functions
c                        of the third kind
c               ms23de : integer vector of exponents corresponding to
c                        ams23dc
c               ace    : array of lnum values of the cosine angular
c                        functions for each of the narg angles
c               aced   : array of lnum values of the first derivatives
c                        of the cosine angular functions for each of
c                        the narg angles
c               ase    : array of lnum values of the sine angular
c                        functions for each of the narg angles
c               ased   : array of lnum values of the first derivatives
c                        of the sine angular functions for each of
c                        the narg angles
c
        use param
c
c  real(knd) scalars
        real(knd) api,a01,asubl,b12,bsubl,cm,cepio2,cedpio2,ce0,dec,
     1            eigaval,eigbval,eigavalp,eigbvalp,eiga1,eiga2,eiga3,
     2            eiga4,eiga5,eigb1,eigb2,eigb3,eigb4,eigb5,esa,esb,
     3            gamma,m1botc,m1bots,mc1c,mc1dc,mc2c,mc2dc,mc3c,mc3dc,
     4            ms1c,ms1dc,ms2c,ms2dc,ms3c,ms3dc,m3bot,pi,q,rl,sed0,
     5            sepio2,sedpio2,sgna,sgnb,wroncc,wroncs,wront,xi,xb,
     6            x1,z
c
c  real(knd) vectors with dimension lnum
        real(knd) eiga(lnum),eigb(lnum)
        real(knd) amc1c(lnum),amc1dc(lnum),amc23c(lnum),amc23dc(lnum),
     1            ams1c(lnum),ams1dc(lnum),ams23c(lnum),ams23dc(lnum)
c
c  integer vectors with dimension lnum
        integer   mc1e(lnum),mc1de(lnum),mc23e(lnum),mc23de(lnum),
     1            ms1e(lnum),ms1de(lnum),ms23e(lnum),ms23de(lnum)
c
c  real(knd) arrays with dimensions lnum and narg
        real(knd) ace(lnum,narg),aced(lnum,narg),ase(lnum,narg),
     1            ased(lnum,narg)
c  
c  real(knd) vectors with dimension maxd
        real(knd) enra(maxd),blista(maxd),glista(maxd),enrb(maxd),
     1            blistb(maxd),glistb(maxd)
c
c  real(knd) vectors with dimension maxj
        real(knd) cbesf(maxj),cbesf1(maxj),cbesf2(maxj),cbesdf(maxj),
     1            cbesdf1(maxj),cbesdf2(maxj),cbesdr(maxj),
     2            cbesdr1(maxj),cbesdr2(maxj)
c
c  integer and real(knd) vectors with dimension maxlp
        dimension ibese(maxlp),ibese1(maxlp),ibese2(maxlp),
     1            ineue(maxlp),ineue2(maxlp),ineue3(maxlp)
        real(knd) cbesn(maxlp),cbesn1(maxlp),cbesn2(maxlp),cneun(maxlp),
     1            cneun2(maxlp),cneun3(maxlp)
c       
c  real(knd) vectors with dimension maxn
        real(knd) cneudf2(maxn),cneuf2(maxn),cneudr2(maxn)
c
c  real(knd) vectors with dimension maxk
        real(knd) cneudf(maxk),cneuf(maxk),cneudr(maxk)

c  real(knd) vectors with dimension maxkbp
        real(knd) cneudf3(maxkbp),cneuf3(maxkbp),cneudr3(maxkbp)
c
c  real(knd) vectors with dimension given by maxp
        real(knd) cosi(narg,maxp),sine(narg,maxp)
c
c  integer and real(knd) vectors with dimension narg
        dimension nacca(narg),naccc(narg),naccs(narg)
        real(knd) arg(narg),barg(narg),ce(narg),ced(narg),
     1            se(narg),sed(narg)
c
c  real(knd) vectors with dimension ngau
        real(knd) xr(ngau),wr(ngau)
c
        dec=10.0e0_knd**(-ndec-1)
        xi=x1+1.0e0_knd
        pi=acos(-1.0e0_knd)
        api=pi/180.0e0_knd
        gamma=0.577215664901532860606512090082402431042e0_knd
        igauss=0
c
c  begin loops
c
            if(iopang.ne.0) then
              do jarg=1,narg
              barg(jarg)=arg(jarg)*api
              end do
            limsc=2*(lnum+cm+2*ndec+50)
            call sincos (limsc,maxp,narg,isq,barg,sine,cosi)
            end if
c            if(izxi.eq.1) then
c            if(ioprad.ne.0.and.knd.eq.kindd) write(20,10) z
c            if(ioprad.ne.0.and.knd.eq.kindd) write(40,10) z
c            if(ioprad.ne.0.and.knd.eq.kindq) write(20,15) z
c            if(ioprad.ne.0.and.knd.eq.kindq) write(40,15) z
c10          format(1x,' z = ',e24.15)
c15          format(1x,' z = ',e40.31)
c            end if
c            if(izxi.eq.2) then
c            if(ioprad.ne.0.and.knd.eq.kindd) write(20,20) xi
c            if(ioprad.ne.0.and.knd.eq.kindd) write(40,20) xi
c            if(ioprad.ne.0.and.knd.eq.kindq) write(20,25) xi
c            if(ioprad.ne.0.and.knd.eq.kindq) write(40,25) xi
c            end if
c20          format(1x,'xi = ',e24.15)
c25          format(1x,'xi = ',e40.31)
c            if(icq.eq.1) then
c            if(ioprad.ne.0.and.knd.eq.kindd) write(20,30) q
c            if(ioprad.ne.0.and.knd.eq.kindd) write(40,30) q
c            if(iopang.ne.0.and.knd.eq.kindd) write(30,30) q
c            if(iopang.ne.0.and.knd.eq.kindd) write(50,30) q
c            if(ioprad.ne.0.and.knd.eq.kindq) write(20,35) q
c            if(ioprad.ne.0.and.knd.eq.kindq) write(40,35) q
c            if(iopang.ne.0.and.knd.eq.kindq) write(30,35) q
c            if(iopang.ne.0.and.knd.eq.kindq) write(50,35) q
c30          format(1x,' q = ',e24.15)
c35          format(1x,' q = ',e40.31)
c            end if
c            if(icq.eq.2) then
c            if(ioprad.ne.0.and.isq.eq.1.and.knd.eq.kindd)
c     1         write(20,40) cm
c            if(ioprad.ne.0.and.isq.eq.1.and.knd.eq.kindd)
c     1         write(40,40) cm
c            if(ioprad.ne.0.and.isq.eq.-1.and.knd.eq.kindd)
c     1         write(20,50) cm
c            if(ioprad.ne.0.and.isq.eq.-1.and.knd.eq.kindd)
c     1         write(40,50) cm
c            if(ioprad.ne.0.and.isq.eq.1.and.knd.eq.kindq)
c     1         write(20,45) cm
c            if(ioprad.ne.0.and.isq.eq.1.and.knd.eq.kindq)
c     1         write(40,45) cm
c            if(ioprad.ne.0.and.isq.eq.-1.and.knd.eq.kindq)
c     1         write(20,55) cm
c            if(ioprad.ne.0.and.isq.eq.-1.and.knd.eq.kindq)
c     1         write(40,55) cm
c40          format(1x,' c = ',e24.15)
c45          format(1x,' c = ',e40.31)
c50          format(1x,' c = i times',e24.15)
c55          format(1x,' c = i times',e40.31)
c            if(iopang.ne.0.and.isq.eq.1.and.knd.eq.kindd)
c     1         write(30,40) cm
c            if(iopang.ne.0.and.isq.eq.1.and.knd.eq.kindd)
c     1         write(50,40) cm
c            if(iopang.ne.0.and.isq.eq.-1.and.knd.eq.kindd)
c     1         write(30,50) cm
c            if(iopang.ne.0.and.isq.eq.-1.and.knd.eq.kindd)
c     1         write(50,50) cm
c            if(iopang.ne.0.and.isq.eq.1.and.knd.eq.kindq)
c     1         write(30,45) cm
c            if(iopang.ne.0.and.isq.eq.1.and.knd.eq.kindq)
c     1         write(50,45) cm
c            if(iopang.ne.0.and.isq.eq.-1.and.knd.eq.kindq)
c     1         write(30,55) cm
c            if(iopang.ne.0.and.isq.eq.-1.and.knd.eq.kindq)
c     1         write(50,55) cm
c            end if
          wront=2.0e0_knd/pi
          ibflag=1
          iopbpe=1
          iopneu=0
          iopbp1=0
          iopbes=1
          if(isq.eq.1) go to 60
          iopneu=1
          if(minacc.le.(ndec-4).and.cm.le.5.0e0_knd) iopneu=0
          if(minacc.le.(ndec-6).and.cm.le.10.0e0_knd) iopneu=0
          if(minacc.le.(ndec-8).and.cm.le.15.0e0_knd) iopneu=0
          if(minacc.le.(ndec-10).and.cm.le.20.0e0_knd) iopneu=0
          if(iopneu.eq.0) go to 60
          iopbpe=0
          iopbes=0
          iopbp1=1
60        if(x1.eq.0.0e0_knd) iopneu=0
          if(x1.eq.0.0e0_knd) iopbpe=1
          if(ioprad.eq.0) go to 70
          limj=lnum+4*ndec+int(cm)+100
          xb=sqrt(x1*(x1+2.0e0_knd))
          limbes=2*int(cm*xb)+2*ndec
          limjt=max(limbes,limj)
          numlp=lnum
          iopr=2
          if(isq.eq.1.and.x1.ne.0.0e0_knd) call cbesj(cm,xb,limj,maxj,
     1                                 numlp,maxlp,limjt,ndec,limbes,
     2                                 cbesf,cbesn,ibese,cbesdf,
     3                                 cbesdr,iopr)
          if(isq.eq.-1.and.x1.ne.0.0e0_knd) call cbesi(cm,xb,limj,maxj,
     1                                 numlp,maxlp,limjt,ndec,limbes,
     2                                 cbesf,cbesn,ibese,cbesdf,
     3                                 cbesdr,iopr,nex)
          xb=0.5e0_knd/(xi+sqrt(x1*(x1+2.0e0_knd)))
          limbes=2*int(cm*xb)+2*ndec
          limbpe=1.5*lnum+4*ndec+int(cm)+100
          numlp=lnum+1
          if(isq.eq.1) go to 65
          ismax=int(0.85e0_knd*cm)
          if(x1.lt.1.0e0_knd) ismax=ismax*1.5
          if(x1.lt.0.1e0_knd.and.x1.ge.0.00001e0_knd) ismax=ismax*
     1                           int(3.0e0_knd**(-1.0e0_knd-log10(x1)))
          if(x1.lt.0.00001e0_knd) ismax=ismax*81
          limbpe=max(limbpe,ismax+int(cm/pi)+4*ndec+
     1                              int(cm)+102)
          numlp=max(numlp,ismax+int(cm/pi)+2)
65        limjt=max(limbes,limbpe)
          iopr=1
          if(isq.eq.1) call cbesj(cm,xb,limbpe,maxj,numlp,maxlp,limjt,
     1                            ndec,limbes,cbesf1,cbesn1,
     2                            ibese1,cbesdf1,cbesdr1,iopr)
          if(isq.eq.-1) call cbesi(cm,xb,limbpe,maxj,numlp,maxlp,limjt,
     1                             ndec,limbes,cbesf1,cbesn1,
     2                             ibese1,cbesdf1,cbesdr1,iopr,nex)
70        eiga1=0.0e0_knd
          eiga2=0.0e0_knd
          eiga3=0.0e0_knd
          eiga4=0.0e0_knd
          eiga5=0.0e0_knd
          eigb1=0.0e0_knd
          eigb2=0.0e0_knd
          eigb3=0.0e0_knd
          eigb4=0.0e0_knd
          eigb5=0.0e0_knd
          kflag=0
          kflaga1=0
          kflagb1=0
          m1flag=0
          jbesa=0
          jbesb=0
          jbesmax=0
          jbp1a=0
          jbp1b=0
          jbp1=0
          jneua=0
          jneub=0
          jbpe=0
          jbpemax=0
          jneumax=0
          iss1=0
          iss2=0
          iflag=0
          eigval=0.0q0
          jmfb=0
          jsubb=0
          jbpeb=0
c          if(ioprad.ne.0.and.isq.eq.1.and.x1.eq.0.0q0) write(40,75)
c75        format(1x,/,5x,'Note that when x = 1.0, some of the radial',
c     1           ' functions of the second kind can be very small for',
c     2           /,5x,'l less than 2c/pi. These are correspondingly',
c     3           ' inaccurate. This should not reduce the accuracy',
c     4           /,5x,'of solutions to problems where the radial',
c     5           ' functions of the second kind appear in a pair with',
c     6           /,5x,'the radial functions of the first kind. The',
c     7           ' radial functions of the first kind are very',
c     8           /,5x,'accurate and dominate the radial functions',
c     9           ' of the second kind to the degree the radial',
c     A           /,5x,'functions of the second kind are inaccurate.',
c     B           ' Note that the Wronskian accuracy overestimates',
c     C           /,5x,'the accuracy of the smaller function values',/)
c          if(ioprad.ne.0.and.isq.eq.-1.and.x1.eq.0.0q0) write(40,80)
c80        format(1x,/,5x,'When x = 1.0 and c is not small and l is',
c     1           ' less than about 2c/pi, accuracy of the radial',
c     2           /,5x,'functions of the third kind calculated using',
c     3           ' the Bessel product expansion is reduced due to',
c     4           /,5x,'subtraction errors. These errors are similar',
c     5           ' in size to those obtained in computing the angular',
c     6           /,5x,'functions ce(0) and sed(0). The associated',
c     7           ' radial functions of the first kind are computed',
c     8           /,5x,'accurately. When the Wronskian accuracy is less',
c     9           ' than the desired minimum accuracy, we use the',
c     A           /,5x,'Wronskian to compute an accurate value for the',
c     B           ' radial function of the third kind (for sine',
c     C           /,5x,'functions). We also do this for the first',
c     D           ' derivative of the radial function of the third kind',
c     E           /,5x,'(for cosine functions). When the Wronskian',
c     F           ' accuracy equals zero digits, the radial functions',
c     G           /,5x,'of the third kind (for cosine functions) and',
c     H           ' their first derivatives (for sine functions) are',
c     I           /,5x,'set  equal to zero.',/)
            do 700 li=1,lnum
            l=li-1
c            if(iopang.ne.0) write(30,90) l
c90          format(1x,'l = ',i6)
            ix=l-2*(l/2)
            rl=real(l,knd)
            istop=0
            naccrs=0
            naccrc=0
            issd=max(iss1,iss2)
            limdrad=4*ndec+int(cm)+100
            if(iopbes.ne.0.and.l.gt.1) limdrad=jbesmax+jbesmax+20+cm/25
            if(iopbpe.ne.0.and.l.gt.1) limdrad=max(limdrad,jbpemax+
     1                                          jbpemax+40+int(cm/25))
            if(iopbp1.ne.0.and.l.gt.1) limdrad=max(limdrad,jbp1+jbp1+
     1                20+int(cm/25)+2*issd,2*l+20+int(cm/25)+2*issd)
            limdang=4*ndec+2*int(cm)+100
            if(iopang.ne.0.and.l.gt.1) limdang=jang+jang+20+cm/25
            if(iopang.eq.0) limd=limdrad
            if(ioprad.eq.0) limd=limdang
            if(iopang.ne.0.and.ioprad.ne.0) limd=max(limdang,limdrad)
            if(l.eq.0) limmf=limdang
            if(l.gt.0) limmf=jmf+jmf+40+cm/25
            limd=max(limd,limmf)
            if(isq.eq.1.or.ioprad.eq.0.or.iopneu.eq.0) go to 100
            if(iopneu.eq.3.and.l.gt.1) go to 95
            limdneu=1.5*cm+100
            if(x1.lt.1.0e0_knd.and.x1.ge.0.01e0_knd) limdneu=40/x1+
     1                                   cm/sqrt(x1)+1.5*cm+200
            if(x1.lt.0.01e0_knd.and.x1.gt.0.0e0_knd) limdneu=34/x1+
     1                                  1.5*cm+1.4*cm/sqrt(x1)+200
95          continue
            if(iopneu.eq.3.and.l.gt.1) limdneu=jneumax+jneumax+limndel+
     1                                         100
            limdneu=min(limdneu,maxk)
            if(iopneu.ne.0) limd=max(limd,limdneu)
100         if(limd.gt.2*maxd) limd=2*maxd-5
            if(l.eq.0) ienra=(3*ndec+int(cm))/2
            call geteiga(l,cm,eiga1,eiga2,eiga3,eiga4,eiga5,eigaval)
            call convera(l,cm,limd,eiga1,eiga3,eiga4,eiga5,ndec,maxd,
     1                    enra,ienra,sgna,a01,ia01,blista,glista,
     2                    eigaval)
            if(l.eq.1) ienrb=(3*ndec+int(cm))/2
            if(l.ne.0) call geteigb(l,cm,eigb1,eigb2,eigb3,eigb4,eigb5,
     1                              eigbval)
            if(l.ne.0) call converb(l,cm,limd,eigb1,eigb3,eigb4,eigb5,
     1                              ndec,maxd,enrb,ienrb,sgnb,b12,ib12,
     2                              blistb,glistb,eigbval)
            esa=eigaval
            if(l.ne.0) esb=eigbval
            if(isq.eq.-1.and.ix.eq.1) esa=eigbval
            if(isq.eq.-1.and.ix.eq.1) esb=eigaval
            eiga(l+1)=esa
            if(l.ne.0) eigb(l+1)=esb
            if(l.eq.0) eigb(l+1)=0.0e0_knd
c              if(knd.eq.kindd) then
c              if(ioprad.ne.0.and.l.eq.0) write(40,110) l,esa
c110           format(1x,'l =',i5,5x,'eigenvalue a =',e24.15)
c              if(ioprad.ne.0.and.l.ne.0) write(40,120) l,esa,esb
c120           format(1x,'l =',i5,5x,'eigenvalue a =',e24.15,/,14x,
c     1               'eigenvalue b =',e24.15)
c              if(iopang.ne.0.and.l.ne.0) write(50,120) l,esa,esb
c              if(iopang.ne.0.and.l.eq.0) write(50,110) l,esa
c              end if
c              if(knd.eq.kindq) then
c              if(ioprad.ne.0.and.l.eq.0) write(40,115) l,esa
c115           format(1x,'l =',i5,5x,'eigenvalue a =',e40.31)
c              if(ioprad.ne.0.and.l.ne.0) write(40,125) l,esa,esb
c125           format(1x,'l =',i5,5x,'eigenvalue a =',e40.31,/,14x,
c     1               'eigenvalue b =',e40.31)
c              if(iopang.ne.0.and.l.ne.0) write(50,125) l,esa,esb
c              if(iopang.ne.0.and.l.eq.0) write(50,115) l,esa
c              end if
            if(ioprad.eq.0) go to 610
            call dnorma(l,limd,maxd,ndec,enra,ce0,
     1                  cepio2,cedpio2,jmfa,jsuba)
            if(l.ne.0) call dnormb(l,limd,maxd,ndec,enrb,sed0,
     1                             sepio2,sedpio2,jmfb,jsubb)
            jmf=max(jmfa,jmfb)
            jsubmax=max(jsuba,jsubb)
            jtest=ndec-jsubmax
            if(jtest.lt.0) jtest=0
            if(isq.eq.1.or.kflag.eq.1) go to 140
            if(jtest.lt.minacc.and.jsubmax.gt.1) go to 140
            kflag=1
            iopbpe=1
            iopbp1=0
            iopneu=0
            iopbes=1
            if(m1flag.eq.0) m1flag=1
140         if(l.gt.0.and.eigaval.le.eigavalp) go to 710
            if(l.gt.1.and.eigbval.le.eigbvalp) go to 710
            eigavalp=eigaval
            if(l.ne.0) eigbvalp=eigbval
            if(ioprad.eq.0) go to 610
c
c  determine Mathieu radial functions of the first kind
c
c    calculation of m1 using a series of cylindrical Bessel functions
c    of the first kind with argument c*sqrt(xi*xi-1)
            if(iopbes.eq.0) go to 170
c            if(x1.ne.0.0e0_knd) write(40,150)
c150         format(1x,'m1 calculation using series of Bessel functions',
c     1             ' with argument cm*sqrt(xi*xi-1)')
c            if(x1.eq.0.0e0_knd) write(40,155)
c155         format(1x,'m1 equal to nonzero term (if any) in series of',
c     1             ' Bessel functions with argument cm*sqrt(xi*xi-1)')
            if(iopbes.eq.1) limr1=l+4*ndec+int(cm)+100
            if(iopbes.eq.2) limr1=jbesa+jbesa+20+cm/25
            jbesa1=jbesa
            iopcs=1
            if(isq.eq.1.and.ix.eq.0) m1botc=cepio2
            if(isq.eq.1.and.ix.eq.1) m1botc=cedpio2
            if(isq.eq.-1.and.ix.eq.0) m1botc=ce0
            if(isq.eq.-1.and.ix.eq.1) m1botc=sed0
            call m1bes(iopcs,l,cm,x1,isq,limr1,ndec,maxd,enra,enrb,
     1                 maxj,maxlp,nex,cbesf,cbesn,ibese,cbesdf,cbesdr,
     2                 m1botc,a01,ia01,b12,ib12,jbesa,kflaga1,nsubc1,
     3                 mc1c,imc1e,mc1dc,imc1de)
            jbesa=max(jbesa1,jbesa)
            if(x1.eq.0.0e0_knd) jbesa=jmfa
            if(l.eq.0) go to 160
            if(iopbes.eq.1.or.l.eq.1) limr1=l+4*ndec+int(cm)+100
            if(iopbes.eq.2.and.l.ne.1) limr1=jbesb+jbesb+20+cm/25
            iopcs=2
            jbesb1=jbesb
            if(isq.eq.1.and.ix.eq.0) m1bots=sedpio2
            if(isq.eq.1.and.ix.eq.1) m1bots=sepio2
            if(isq.eq.-1.and.ix.eq.0) m1bots=sed0
            if(isq.eq.-1.and.ix.eq.1) m1bots=ce0
            call m1bes(iopcs,l,cm,x1,isq,limr1,ndec,maxd,enra,enrb,
     1                 maxj,maxlp,nex,cbesf,cbesn,ibese,cbesdf,cbesdr,
     2                 m1bots,a01,ia01,b12,ib12,jbesb,kflagb1,nsubs1,
     3                 ms1c,ims1e,ms1dc,ims1de)
            jbesb=max(jbesb1,jbesb)
            if(x1.eq.0.0e0_knd) jbesb=jmfb
160         iopbes=2
            jbesmax=max(jbesa,jbesb)
            if(isq.eq.1) go to 200
c
c    calculation of m1 using a series of products of cylindrical
c    Bessel functions of the first kind (q < 0, c imaginary)
170         if(iopbp1.eq.0) go to 200
c            write(40,180)
c180         format(1x,'m1 calculation using Bessel product expansion')
            if(l.ne.0.or.x1.eq.0.0e0_knd) go to 190
            limbpe=lnum+4*ndec+int(cm)+100
            numi=max(limbpe,ismax+int(cm/pi)+4*ndec+int(cm)+102)
            numlp=ismax+int(cm/pi)+2
            xb=0.5e0_knd*(xi+sqrt(x1*(x1+2.0e0_knd)))
            limbes=2*int(cm*xb)+2*ndec
            limjt=max(limbes,numi)
            iopr=1
            call cbesi(cm,xb,numi,maxj,numlp,maxlp,limjt,ndec,limbes,
     1                 cbesf2,cbesn2,ibese2,cbesdf2,cbesdr2,iopr,nex)
190         continue
            jbp1a1=jbp1a
            iopcs=1
            if(l.eq.0) limbp1a=4*ndec+int(cm)+100
            if(l.gt.0) limbp1a=jbp1+jbp1+20+cm/25
            if(x1.eq.0.0e0_knd) call m1bpe0(iopcs,iss1,ismax,l,cm,
     1                               limbp1a,ndec,limd,maxd,enra,enrb,
     2                               maxj,maxlp,nex,cbesf1,cbesn1,
     3                               ibese1,cbesdf1,cbesdr1,a01,ia01,
     4                               b12,ib12,jbp1a,nsubc1,mc1c,imc1e,
     5                               mc1dc,imc1de)
            if(x1.ne.0.0e0_knd) call m1bpe(iopcs,iss1,ismax,l,cm,x1,
     1                               limbp1a,ndec,limd,maxd,enra,enrb,
     2                               maxj,maxlp,nex,cbesf1,cbesn1,
     3                               ibese1,cbesdf1,cbesdr1,cbesf2,
     4                               cbesn2,ibese2,cbesdf2,cbesdr2,a01,
     5                               ia01,b12,ib12,jbp1a,nsubc1,mc1c,
     6                               imc1e,mc1dc,imc1de)
            jbp1a=max(jbp1a1,jbp1a)
            if(l.eq.0) go to 200
            iopcs=2
            jbp1b1=jbp1b
            if(l.eq.1) limbp1a=4*ndec+int(cm)+100
            if(l.gt.1) limbp1a=jbp1+jbp1+20+cm/25
            if(x1.eq.0.0e0_knd) call m1bpe0(iopcs,iss2,ismax,l,cm,
     1                               limbp1a,ndec,limd,maxd,enra,enrb,
     2                               maxj,maxlp,nex,cbesf1,cbesn1,
     3                               ibese1,cbesdf1,cbesdr1,a01,ia01,
     4                               b12,ib12,jbp1b,nsubs1,ms1c,ims1e,
     5                               ms1dc,ims1de)
            if(x1.ne.0.0e0_knd) call m1bpe(iopcs,iss2,ismax,l,cm,x1,
     1                               limbp1a,ndec,limd,maxd,enra,enrb,
     2                               maxj,maxlp,nex,cbesf1,cbesn1,
     3                               ibese1,cbesdf1,cbesdr1,cbesf2,
     4                               cbesn2,ibese2,cbesdf2,cbesdr2,a01,
     5                               ia01,b12,ib12,jbp1b,nsubs1,ms1c,
     6                               ims1e,ms1dc,ims1de)
            jbp1b=max(jbp1b1,jbp1b)
200         jbp1=max(jbp1a,jbp1b)
c              if(knd.eq.kindd) then
c              if(isq.eq.1.or.(isq.eq.-1.and.ix.eq.0))
c     1             write(40,210) mc1c,imc1e,mc1dc,imc1de
c210           format(10x,'mc1 = ',f17.14,i8,5x,'mc1d = ',f17.14,i8)
c              if((isq.eq.1.or.(isq.eq.-1.and.ix.eq.0)).and.l.ne.0)
c     1             write(40,220) ms1c,ims1e,ms1dc,ims1de
c220           format(10x,'ms1 = ',f17.14,i8,5x,'ms1d = ',f17.14,i8)
c              if(isq.eq.-1.and.ix.eq.1)
c     1             write(40,230) mc1c,imc1e,mc1dc,imc1de
c230           format(10x,'mc1 = i times ',f17.14,i8,5x,
c     1               'mc1d = i times ',f17.14,i8)
c              if(isq.eq.-1.and.ix.eq.1.and.l.ne.0)
c     1             write(40,240) ms1c,ims1e,ms1dc,ims1de
c240           format(10x,'ms1 = i times ',f17.14,i8,5x,
c     1               'ms1d = i times ',f17.14,i8)
c              end if
c              if(knd.eq.kindq) then
c              if(isq.eq.1.or.(isq.eq.-1.and.ix.eq.0))
c     1             write(40,215) mc1c,imc1e,mc1dc,imc1de
c215           format(10x,'mc1 = ', f33.30,i8,5x,'mc1d = ',f33.30,i8)
c              if((isq.eq.1.or.(isq.eq.-1.and.ix.eq.0)).and.l.ne.0)
c     1             write(40,225) ms1c,ims1e,ms1dc,ims1de
c225           format(10x,'ms1 = ', f33.30,i8,5x,'ms1d = ',f33.30,i8)
c              if(isq.eq.-1.and.ix.eq.1)
c     1             write(40,235) mc1c,imc1e,mc1dc,imc1de
c235           format(10x,'mc1 = i times ',f33.30,i8,5x,
c     1               'mc1d = i times ',f33.30,i8)
c              if(isq.eq.-1.and.ix.eq.1.and.l.ne.0)
c     1             write(40,245) ms1c,ims1e,ms1dc,ims1de
c245           format(10x,'ms1 = i times ',f33.30,i8,5x,
c     1               'ms1d = i times ',f33.30,i8)
c              end if
            continue
c
c  determine Mathieu radial functions of the second kind (when q
c  is positive or Mathieu radial functions of the third kind when
c  q is negative)
c
c    calculation of m2 (q positive) or m3 (q negative) using a series
c    of products of cylindrical Bessel functions of the first and second
c    kinds
c
            if(iopbpe.eq.0) go to 400
c            if(isq.eq.1) write(40,250)
c250         format(1x,'m2 calculation using Bessel product series')
c            if(isq.eq.-1) write(40,260)
c260         format(1x,'m3 calculation using Bessel product series')
            limbpe=lnum+4*ndec+int(cm)+100
            if(iopbpe.gt.1) go to 280
            xb=0.5e0_knd*(xi+sqrt(x1*(x1+2.0e0_knd)))
            iopr=1
            limbes=2*int(cm*xb)+2*ndec
            numlp=lnum
            if(igauss.eq.1.or.isq.ne.-1) go to 270
            if((cm*xb.le.1.0e0_knd).or.(cm*xb.ge.40.0e0_knd)) go to 270
            call gauss (ngau,ndec,xr,wr)
            igauss=1
270         if(isq.eq.1) call cneuy (cm,xb,limbpe,maxn,numlp,maxlp,ndec,
     1                               limbes,pi,gamma,cneuf2,cneun2,
     2                               ineue2,cneudf2,cneudr2,iopr)
            if(isq.eq.-1) call cneuk (cm,xb,limbpe,maxkbp,numlp,maxlp,
     1                                ndec,nex,limbes,pi,gamma,ngau,wr,
     2                                xr,iopr,cneuf3,cneun3,ineue3,
     3                                cneudf3,cneudr3)
            iopbpe=iopbpe+1
280         continue
            if(iopbpe.eq.2) limbpea=l+4*ndec+int(cm)+100
            if(iopbpe.eq.3) limbpea=jbpea+jbpea+40+cm/25
            if(limbpea.gt.limbpe) limbpea=limbpe
            iopcs=1
            if(isq.eq.-1) go to 290
            call m2bpe (iopcs,iss1,l,cm,x1,limbpea,ndec,maxd,enra,enrb,
     1                  maxj,maxlp,maxn,nex,cbesf1,cbesn1,ibese1,
     2                  cbesdf1,cbesdr1,cneuf2,cneun2,ineue2,cneudf2,
     3                  cneudr2,a01,ia01,jbpea,nsubbpec,mc2c,imc2e,
     4                  mc2dc,imc2de)
            wroncc=mc1c*mc2dc*10.0e0_knd**(imc1e+imc2de)
            if(x1.ne.0.0e0_knd) wroncc=wroncc-
     1                             mc2c*mc1dc*10.0e0_knd**(imc2e+imc1de)
            naccrc=-int(log10(abs((wroncc-(2.0e0_knd/pi))/
     1              (2.0e0_knd/pi))+dec))
            if(naccrc.lt.0) naccrc=0
            naccrc=min(naccrc,ndec-1)
            go to 300
290         call m3bpe (iopcs,l,cm,x1,limbpea,ndec,maxd,enra,enrb,
     1                  maxj,maxlp,maxkbp,nex,cbesf1,cbesn1,ibese1,
     2                  cbesdf1,cbesdr1,cneuf3,cneun3,ineue3,cneudf3,
     3                  cneudr3,a01,ia01,b12,ib12,pi,jbpea,nsubc3,
     4                  mc3c,imc3e,mc3dc,imc3de)
            wroncc=mc1c*mc3dc*10.0e0_knd**(imc1e+imc3de)
            if(x1.ne.0.0e0_knd) wroncc=wroncc-
     1                             mc3c*mc1dc*10.0e0_knd**(imc3e+imc1de)
            naccrc=-int(log10(abs((wroncc-(2.0e0_knd/pi))/
     1              (2.0e0_knd/pi))+dec))
            if(naccrc.lt.0) naccrc=0
            naccrc=min(naccrc,ndec-1)
            if(x1.ne.0.0e0_knd) go to 300
            if(naccrc.eq.0) mc3c=0
            if(naccrc.eq.0) imc3e=0
            if(naccrc.lt.minacc) mc3dc=(2.0e0_knd/pi)/mc1c
            iterm=int(log10(abs(mc3dc)))
            imc3de=-imc1e+iterm
            mc3dc=mc3dc*10.0e0_knd**(-iterm)
            if(abs(mc3dc).ge.1.0e0_knd) go to 300
            mc3dc=mc3dc*10.0e0_knd
            imc3de=imc3de-1
300         if(l.eq.0) go to 320
            limbpeb=limbpea
            if(l.gt.1.and.iopbpe.eq.3) limbpeb=jbpeb+jbpeb+40+cm/25
            if(limbpeb.gt.limbpe) limbpeb=limbpe
            iopcs=2
            if(isq.eq.-1) go to 310
            call m2bpe (iopcs,iss2,l,cm,x1,limbpeb,ndec,maxd,enra,enrb,
     1                  maxj,maxlp,maxn,nex,cbesf1,cbesn1,ibese1,
     2                  cbesdf1,cbesdr1,cneuf2,cneun2,ineue2,cneudf2,
     3                  cneudr2,b12,ib12,jbpeb,nsubbpes,ms2c,ims2e,
     4                  ms2dc,ims2de)
            wroncs=-ms2c*ms1dc*10.0e0_knd**(ims2e+ims1de)
            if(x1.ne.0.0e0_knd) wroncs=wroncs+
     1                             ms1c*ms2dc*10.0e0_knd**(ims1e+ims2de)
            naccrs=-int(log10(abs((wroncs-(2.0e0_knd/pi))/
     1              (2.0e0_knd/pi))+dec))
            if(naccrs.lt.0) naccrs=0
            naccrs=min(naccrs,ndec-1)
            go to 320
310         call m3bpe (iopcs,l,cm,x1,limbpeb,ndec,maxd,enra,enrb,
     1                  maxj,maxlp,maxkbp,nex,cbesf1,cbesn1,ibese1,
     2                  cbesdf1,cbesdr1,cneuf3,cneun3,ineue3,cneudf3,
     3                  cneudr3,a01,ia01,b12,ib12,pi,jbpeb,nsubs4,
     4                  ms3c,ims3e,ms3dc,ims3de)   
            wroncs=-ms3c*ms1dc*10.0e0_knd**(ims3e+ims1de)
            if(x1.ne.0.0e0_knd) wroncs=wroncs+
     1                             ms1c*ms3dc*10.0e0_knd**(ims1e+ims3de)
            naccrs=-int(log10(abs((wroncs-(2.0e0_knd/pi))/
     1              (2.0e0_knd/pi))+dec))
            if(naccrs.lt.0) naccrs=0
            naccrs=min(naccrs,ndec-1)
            if(x1.ne.0.0e0_knd) go to 320
            if(naccrs.eq.0) ms3dc=0.0e0_knd
            if(naccrs.eq.0) ims3de=0
            if(naccrs.lt.minacc) ms3c=-(2.0e0_knd/pi)/ms1dc
            iterm=int(log10(abs(ms3c)))
            ims3e=-ims1de+iterm
            ms3c=ms3c*10.0e0_knd**(-iterm)
            if(abs(ms3c).ge.1.0e0_knd) go to 320
            ms3c=ms3c*10.0e0_knd
            ims3e=ims3e-1
320         if(l.eq.0) jbpeb=0
            jbpemax=max(jbpea,jbpeb)
            iopbpe=3
            if(isq.eq.-1) go to 350
c              if(knd.eq.kindd) then
c              write(40,330) naccrc,mc2c,imc2e,mc2dc,imc2de
c330           format(15x,'Wronskian accuracy =',i3,
c     1               ' decimal digits.'/,10x,'mc2 = ', f17.14,i6,
c     2               5x,'mc2d = ',f17.14,i6)
c              if(l.gt.0) write(40,340) naccrs,ms2c,ims2e,ms2dc,ims2de
c340           format(15x,'Wronskian accuracy =',i3,
c     1               ' decimal digits.'/,10x,'ms2 = ', f17.14,i6,
c     2               5x,'ms2d = ',f17.14,i6)
c              go to 400
c              end if
c              if(knd.eq.kindq) then
c              write(40,335) naccrc,mc2c,imc2e,mc2dc,imc2de
c335           format(15x,'Wronskian accuracy =',i3,
c     1               ' decimal digits.'/,10x,'mc2 = ', f33.30,i6,
c     2               5x,'mc2d = ',f33.30,i6)
c              if(l.gt.0) write(40,345) naccrs,ms2c,ims2e,ms2dc,ims2de
c345           format(15x,'Wronskian accuracy =',i3,
c     1               ' decimal digits.'/,10x,'ms2 = ', f33.30,i6,
c     2               5x,'ms2d = ',f33.30,i6)
c              go to 400
c              end if
350         continue
c              if(knd.eq.kindd) then
c              if(ix.eq.1) write(40,360) naccrc,mc3c,imc3e,mc3dc,imc3de
c360           format(15x,'Wronskian accuracy =',i3,
c     1               ' decimal digits.'/,10x,'mc3 = ', f17.14,i8,
c     2               5x,'mc3d = ',f17.14,i8)
c              if(ix.eq.1.and.l.gt.0) write(40,370) naccrs,ms3c,ims3e,
c     1                                             ms3dc,ims3de
c370           format(15x,'Wronskian accuracy =',i3,
c     1               ' decimal digits.'/,10x,'ms3 = ', f17.14,i8,
c     2               5x,'ms3d = ',f17.14,i8)
c              if(ix.eq.0) write(40,380) naccrc,mc3c,imc3e,mc3dc,imc3de
c380           format(15x,'Wronskian accuracy =',i3,
c     1               ' decimal digits.'/,10x,'mc3 = i times ',f17.14,i8,
c     2               5x,'mc3d = i times ',f17.14,i8)
c              if(l.gt.0.and.ix.eq.0) write(40,390) naccrs,ms3c,ims3e,
c     1                                             ms3dc,ims3de
c390           format(15x,'Wronskian accuracy =',i3,
c     1               ' decimal digits.'/,10x,'ms3 = i times ',f17.14,i8,
c     2               5x,'ms3d = i times ',f17.14,i8)
c              end if
c              if(knd.eq.kindq) then
c              if(ix.eq.1) write(40,365) naccrc,mc3c,imc3e,mc3dc,imc3de
c365           format(15x,'Wronskian accuracy =',i3,
c     1               ' decimal digits.'/,10x,'mc3 = ', f33.30,i8,
c     2               5x,'mc3d = ',f33.30,i8)
c              if(ix.eq.1.and.l.gt.0) write(40,375) naccrs,ms3c,ims3e,
c     1                                             ms3dc,ims3de
c375           format(15x,'Wronskian accuracy =',i3,
c     1               ' decimal digits.'/,10x,'ms3 = ', f33.30,i8,
c     2               5x,'ms3d = ',f33.30,i8)
c              if(ix.eq.0) write(40,385) naccrc,mc3c,imc3e,mc3dc,imc3de
c385           format(15x,'Wronskian accuracy =',i3,
c     1               ' decimal digits.'/,10x,'mc3 = i times ',f33.30,i8,
c     2               5x,'mc3d = i times ',f33.30,i8)
c              if(l.gt.0.and.ix.eq.0) write(40,395) naccrs,ms3c,ims3e,
c     1                                             ms3dc,ims3de
c395           format(15x,'Wronskian accuracy =',i3,
c     1               ' decimal digits.'/,10x,'ms3 = i times ',f33.30,i8,
c     2               5x,'ms3d = i times ',f33.30,i8)
c              end if
400         continue
            if(isq.eq.1) go to 510
c
c    calculation of m3 (q negative) using the conventional Neumann
c    expansion in a series of modified cylindrical Bessel functions
c    of the second kind with argument cm*xi
c
            if(iopneu.eq.0) go to 510
            if(ibflag.eq.0) go to 420
            limn=1.5*cm+100
            if(x1.lt.1.0e0_knd.and.x1.ge.0.01e0_knd) limn=40/x1+cm/
     1                                           sqrt(x1)+1.5*cm+200
            if(x1.lt.0.01e0_knd.and.x1.gt.0.0e0_knd) limn=34/x1+1.5*cm+
     1                                           1.4*cm/sqrt(x1)+200
            limn=min(limn,maxk-3)
            limndel=2*limn/cm
            xb=xi
            limbes=2*int(cm*xb)+2*ndec
            iopr=2
            numlp=min(lnum,int(2.0e0_knd*cm/pi)+10)
            if(igauss.eq.1) go to 410
            if((cm*xb.le.1.0e0_knd).or.(cm*xb.ge.40.0e0_knd)) go to 410
            call gauss (ngau,ndec,xr,wr)
            igauss=1
410         call cneuk(cm,xb,limn,maxk,numlp,maxlp,ndec,nex,limbes,pi,
     1                 gamma,ngau,wr,xr,iopr,cneuf,cneun,ineue,cneudf,
     2                 cneudr)
            ibflag=0
420         if(iopneu.eq.1) iopneu=2
            if(iopneu.eq.4) go to 510
            if(iopneu.eq.3.and.l.gt.1) go to 430
            limneua=1.5*cm+100
            if(x1.lt.1.0e0_knd.and.x1.ge.0.01e0_knd) limneua=40/x1+cm/
     1                                      sqrt(x1)+1.5*cm+200
            if(x1.lt.0.01e0_knd.and.x1.gt.0.0e0_knd) limneua=34/x1+
     1                                      1.5*cm+1.4*cm/sqrt(x1)+200
430         continue
            if(iopneu.eq.3.and.l.gt.1) limneua=jneumax+jneumax+limndel+
     1                                 100
            limneua=min(limneua,limn-2)
            iopcs=1
            if(ix.eq.0) m3bot=cepio2
            if(ix.eq.1) m3bot=sepio2
c            write(40,440)
c440         format(1x,'m3 calculation using Bessel function series')
            call m3neu (iopcs,l,cm,x1,limneua,ndec,minacc,maxd,
     1                  enra,enrb,maxk,maxlp,nex,cneuf,cneun,ineue,
     2                  cneudf,cneudr,m3bot,pi,jneua,ndig,mc3c,imc3e,
     3                  mc3dc,imc3de,jtermflag)
            wroncc=mc1c*mc3dc*10.0e0_knd**(imc1e+imc3de)-
     1                    mc3c*mc1dc*10.0e0_knd**(imc3e+imc1de)
            naccrc=-int(log10(abs((wroncc-(2.0e0_knd/pi))/
     1              (2.0e0_knd/pi))+dec))
            if(naccrc.lt.0) naccrc=0
            naccrc=min(naccrc,ndec-1)
            if(l.eq.0) go to 450
            if(jtermflag.eq.1) jneua=jneua+2*limndel
            limneub=limneua
            if(jtermflag.eq.1) limneub=limneub+2*limndel
            if(l.gt.1.and.iopneu.eq.3) limneub=jneumax+jneumax+limndel+
     1                                         100
            limneub=min(limneub,limn-2)
            iopcs=2
            if(ix.eq.0) m3bot=sedpio2
            if(ix.eq.1) m3bot=cedpio2
            if(l.gt.0) call m3neu (iopcs,l,cm,x1,limneub,ndec,
     1                             minacc,maxd,enra,enrb,maxk,maxlp,
     2                             nex,cneuf,cneun,ineue,cneudf,
     3                             cneudr,m3bot,pi,jneub,ndig,ms3c,
     4                             ims3e,ms3dc,ims3de,jtermflag)
            wroncs=ms1c*ms3dc*10.0e0_knd**(ims1e+ims3de)-ms3c*ms1dc*
     1             10.0e0_knd**(ims3e+ims1de)
            naccrs=-int(log10(abs((wroncs-(2.0e0_knd/pi))/
     1              (2.0e0_knd/pi))+dec))
            if(naccrs.lt.0) naccrs=0
            naccrs=min(naccrs,ndec-1)
            if(jtermflag.eq.1) jneub=jneub+2*limndel
450         jneumax=max(jneua,jneub,jneumax,limn-2)
            iopneu=3
c            write(40,460) naccrc
c460         format(15x,'Wronskian accuracy =',i3,' decimal digits.')
c              if(knd.eq.kindd) then
c              if(ix.eq.0) write(40,470) mc3c,imc3e,mc3dc,imc3de
c470           format(10x,'mc3 = i times ',f17.14,i8,
c     1               5x,'mc3d = i times ',f17.14,i8)
c              if(ix.eq.1) write(40,480) mc3c,imc3e,mc3dc,imc3de
c480           format(10x,'mc3 = ', f17.14,i8,5x,'mc3d = ',f17.14,i8)
c              if(l.gt.0) write(40,460) naccrs
c              if(ix.eq.0.and.l.gt.0) write(40,490) ms3c,ims3e,ms3dc,
c     1                                             ims3de
c490           format(10x,'ms3 = i times ',f17.14,i8,
c     1               5x,'ms3d = i times ',f17.14,i8)
c              if(ix.eq.1.and.l.gt.0) write(40,500) ms3c,ims3e,ms3dc,
c     1                                             ims3de
c500           format(10x,'ms3 = ', f17.14,i8,5x,'ms3d = ',f17.14,i8)
c              end if
c              if(knd.eq.kindq) then
c              if(ix.eq.0) write(40,475) mc3c,imc3e,mc3dc,imc3de
c475           format(10x,'mc3 = i times ',f33.30,i8,
c     1               5x,'mc3d = i times ',f33.30,i8)
c              if(ix.eq.1) write(40,485) mc3c,imc3e,mc3dc,imc3de
c485           format(10x,'mc3 = ', f33.30,i8,5x,'mc3d = ',f33.30,i8)
c              if(l.gt.0) write(40,460) naccrs
c              if(ix.eq.0.and.l.gt.0) write(40,495) ms3c,ims3e,ms3dc,
c     1                                             ims3de
c495           format(10x,'ms3 = i times ',f33.30,i8,
c     1               5x,'ms3d = i times ',f33.30,i8)
c              if(ix.eq.1.and.l.gt.0) write(40,505) ms3c,ims3e,ms3dc,
c     1                                             ims3de
c505           format(10x,'ms3 = ', f33.30,i8,5x,'ms3d = ',f33.30,i8)
c              end if
510         if(iopneu.eq.4) iopneu=1
            if(isq.eq.-1) go to 540
c            write(20,520) l,mc1c,imc1e,mc1dc,imc1de,mc2c,imc2e,mc2dc,
c     1                    imc2de,naccrc
c520         format(1x,i5,2x,4(f17.14,i6,2x),i2)
c            if(l.gt.0) write(20,530) ms1c,ims1e,ms1dc,ims1de,ms2c,
c     1                               ims2e,ms2dc,ims2de,naccrs
c530         format(8x,4(f17.14,i6,2x),i2)
            go to 600
540         continue            
c            write(20,550) l,mc1c,imc1e,mc1dc,imc1de,mc3c,imc3e,mc3dc,
c     1                    imc3de,naccrc
c550         format(1x,i5,2x,4(f17.14,i7,2x),i2)
c            if(l.gt.0) write(20,560) ms1c,ims1e,ms1dc,ims1de,ms3c,
c     1                               ims3e,ms3dc,ims3de,naccrs
c560         format(8x,4(f17.14,i7,2x),i2)
600         continue
              if(ioprad.eq.1.and.naccrc.lt.6) then
              write(60,*) ' est. acc. = ',naccrc, ' digits for xi = ',
     1                    xi,' cm, q = ',cm,q,' l = ',l
              end if
              if(ioprad.eq.1.and.l.ne.0.and.naccrs.lt.6) then
              write(60,*) ' est. acc. = ',naccrs, ' digits for xi = ',
     1                    xi,' cm, q = ',cm,q,' l = ',l
              end if
              amc1c(li)=mc1c
            amc1dc(li)=mc1dc
            mc1e(li)=imc1e
            mc1de(li)=imc1de
              if(isq.eq.1) then
              amc23c(li)=mc2c
              amc23dc(li)=mc2dc
              mc23e(li)=imc2e
              mc23de(li)=imc2de
              else
              amc23c(li)=mc3c
              amc23dc(li)=mc3dc
              mc23e(li)=imc3e
              mc23de(li)=imc3de
              end if
                if(l.ne.0) then
                ams1c(li)=ms1c
                ams1dc(li)=ms1dc
                ms1e(li)=ims1e
                ms1de(li)=ims1de
                if(isq.eq.1) then
                ams23c(li)=ms2c
                ams23dc(li)=ms2dc
                ms23e(li)=ims2e
                ms23de(li)=ims2de
                else
                ams23c(li)=ms3c
                ams23dc(li)=ms3dc
                ms23e(li)=ims3e
                ms23de(li)=ims3de
                end if
              end if
610         if(iopang.eq.0) go to 700
c
c  determine Mathieu angular functions of the first kind
c
            lims1=4+ndec+2*int(cm)+100
            if(l.ne.0) lims1=jang+jang+20+cm/25
            if(lims1.gt.maxp) lims1=maxp
            call cese(cm,isq,l,lims1,ndec,iopang,maxd,maxp,enra,enrb,
     1                sgna,sgnb,narg,cosi,sine,ce,ced,se,sed,naccc,
     2                naccs,jangc,jangs,asubl,bsubl)
            jang=max(jangc,jangs)
              do 680 jarg=1,narg
              nacca(jarg)=min(naccc(jarg),naccs(jarg))
c              if(iopang.eq.1) write(30,620) arg(jarg),
c     1                       ce(jarg),se(jarg),nacca(jarg)
c              if(iopang.eq.2) write(30,630) arg(jarg),
c     1                       ce(jarg),ced(jarg),se(jarg),
c     2                       sed(jarg),nacca(jarg)
                if(iopang.eq.1) then
                ace(li,jarg)=ce(jarg)
                ase(li,jarg)=se(jarg)
                end if
                if(iopang.eq.2) then
                ace(li,jarg)=ce(jarg)
                aced(li,jarg)=ced(jarg)
                ase(li,jarg)=se(jarg)
                ased(li,jarg)=sed(jarg)
                end if
c620           format(1x,f20.14,5x,e24.15,2x,e24.15,2x,i2)
c630           format(1x,f20.14,5x,e24.15,2x,e24.15,2x,/,26x,
c     1               e24.15,2x,e24.15,2x,i2)
c                if(knd.eq.kindd) then
c                if(iopang.eq.1) write(50,640) arg(jarg),
c     1                                      ce(jarg),se(jarg)
c                if(iopang.eq.2) write(50,650) arg(jarg),
c     1                          ce(jarg),ced(jarg),se(jarg),sed(jarg)
c640             format(5x,f20.14,' degrees',/,10x,'ce = ',
c     1                 e24.15,2x,' se = ',e24.15)
c650             format(5x,f20.14,' degrees',/,10x,'ce = ',e24.15,2x,
c     1                 ' ced = ',e24.15,/,10x,'se = ',e24.15,2x,
c     2                 ' sed = ',e24.15)
c                end if
c                if(knd.eq.kindq) then
c                if(iopang.eq.1) write(50,660) arg(jarg),
c     1                     ce(jarg),se(jarg)
c                if(iopang.eq.2) write(50,670) arg(jarg),
c     1                          ce(jarg),ced(jarg),se(jarg),sed(jarg)
c660             format(5x,f20.14,' degrees',/,10x,'ce = ',
c     1                 e40.31,2x,' se = ',e40.31)
c670             format(5x,f20.14,' degrees',/,10x,'ce = ',e40.31,2x,
c     1                 ' ced = ',e40.31,/,10x,'se = ',e40.31,2x,
c     2                 ' sed = ',e40.31)
c                end if
680           continue
700         continue
710       continue
        return
        end
c
c
        subroutine cese (cm,isq,l,lims1,ndec,iopang,maxd,maxp,enra,
     1                   enrb,sgna,sgnb,narg,cosi,sine,ce,ced,se,sed,
     2                   naccc,naccs,jangc,jangs,asubl,bsubl)
c
c  purpose: to calculate the Mathieu angular functions
c           ce and se and their first derivatives with respect to phi
c
c  parameters :
c
c     input   : cm     : magnitude of c
c               isq    : integer = +1 if q is positive (c is real);
c                        = -1 if q is negative (c is imaginary)
C               l      : l
c               lims1  : twice the number of terms available
c               ndec   : number of decimal digits available in
c                        real(knd) arithmetic
c               iopang : integer = 1 if no derivatives of the angular
c                        functions are desired, = 2 otherwise
c               maxd   : size of enr vector
c               maxp   : number of sine and cosine values available
c                        for each desired angle
c               enra   : a coefficient ratios
c               enrb   : b coefficient ratios
c               sgna   : sign of the a coefficient with index = l
c               sgnb   : sign of the b coefficient with index = l
c               narg   : number of angles desired
c               cosi   : array containing maxp cosine values
c                        for each of the narg angles
c               sine   : array containing maxp sine values
c                        for each of the narg angles
c
c     output  : ce     : vector of dimension narg containing values
c                        of the cosine angular functions for each of
c                        the desired angles
c               ced    : corresponding vector for first derivatives
c                        of the cosine angular functions
c               se     : vector of dimension narg containing values
c                        of the sine angular functions for each of
c                        the desired angles
c               sed    : corresponding vector for first derivatives
c                        of the sine angular functions
c               naccs  : vector of dimension narg containing estimates
c                        of the minimum accuracy in decimal digits of
c                        the se and sed function values for each angle
c               naccc  : vector of dimension narg containing estimates
c                        of the minimum accuracy in decimal digits of
c                        the se and sed function values for each angle
c               jangc  : largest number of terms taken for ce
c               jangs  : largest number of terms taken for se
c               asubl  : a coefficient with subscript l
c               bsubl  : b coefficient with subscript l 
c
        use param
c
c  real(knd) scalars and vectors
        real(knd) adec,asubl,bsubl,cm,coef,csum,cdsum,dcon,dec,dnew,
     1            dnewd,dnuma,dnumb,dold,doldd,dterm,fterm,sgna,sgnb,
     2            ssum,sdsum,suma,temp,tempd
        real(knd) ce(narg),ced(narg),cosi(narg,maxp),enra(maxd),
     1            enrb(maxd),se(narg),sed(narg),sine(narg,maxp)
c  integer vectors
        dimension naccc(narg),naccs(narg)
c
        dec=10.0e0_knd**(-ndec-1)
        dcon=dec
        adec=1000.0e0_knd*dec
        csq=cm*cm
        l2=l/2
        lx2=l2-2*(l2/2)
        ix=l-2*l2
        ixx=ix-1
        ixx2=ixx+2
        lim=lims1/2-ix
c
c  compute dnuma, the magnitude of the a coeficient with index = l
        suma=1.0e0_knd
        if(l.eq.0) suma=2.0e0_knd
        coef=1.0e0_knd
        jlow=l+2
        jterm=l2
          do 110 j=jlow,lims1,2
          jterm=jterm+1
          coef=coef*enra(jterm)*enra(jterm)
          suma=suma+coef
          if(abs(coef/suma).lt.dcon) go to 120
110       continue
120     jlow=l
        jn=jterm
        if(jlow.lt.2) go to 140
        coef=1.0e0_knd
        jterm=l2
        j=jlow
          do 130 jj=2,jlow,2
          coef=coef/(enra(jterm)*enra(jterm))
          jterm=jterm-1
          j=j-2
          suma=suma+coef
          if(abs(coef/suma).lt.dcon) go to 140
130       continue
        if(ix.eq.0) suma=suma+coef
140     dnuma=1.0e0_knd/sqrt(suma)
        asubl=dnuma*sgna
c        if(isq.eq.1.or.ix.eq.0) write(50,145) jn,lim
c        if(isq.eq.-1.and.ix.eq.1) write(50,345) jn,lim
c145     format(5x,'ce normalization series converged in ',
c     1         i6,' terms; ',i6,' terms available')
        jlow=l2+1
        jangc=0
          do 290 k=1,narg
c
c  compute the angular function ce
          if(abs(cosi(k,3-ix)).gt.adec) go to 150
          ce(k)=0.0e0_knd
          naccc(k)=ndec
          jfun=1
c          if(iopang.eq.1.and.(isq.eq.1.or.ix.eq.0)) write(50,175) k,jfun
c          if(iopang.eq.1.and.isq.eq.-1.and.ix.eq.1) write(50,375) k,jfun
          go to 210
150       dold=1.0e0_knd
          fterm=0.0e0_knd
          csum=cosi(k,l+1)
          if(csum.gt.0.0e0_knd) fterm=csum
            do 160 j=jlow,lim
            dnew=dold*enra(j)
            dterm=dnew*cosi(k,j+j+ixx2)
            csum=csum+dterm
            if(dterm.gt.0.0e0_knd) fterm=fterm+dterm
            if((abs(dterm/csum).lt.dcon).and.(abs(cosi(k,j+j+ixx2))
     1              .gt.1.0e-01_knd)) go to 170
            dold=dnew
160         continue
170       if(j.gt.jangc) jangc=j
          jfun=j
c          if(iopang.eq.1.and.(isq.eq.1.or.ix.eq.0)) write(50,175) k,jfun
c          if(iopang.eq.1.and.isq.eq.-1.and.ix.eq.1) write(50,375) k,jfun
c175       format(8x,'angle ',i3,'. ce converged in ',i6,' terms')
          if(l2.lt.1) go to 190
          dold=1.0e0_knd
          j=l2
            do 180 jj=1,l2
            dnew=dold/enra(j)
            dterm=dnew*cosi(k,j+j+ixx)
            csum=csum+dterm
            if(dterm.gt.0.0e0_knd) fterm=fterm+dterm
            if(csum.eq.0.0e0_knd) go to 180
            if((abs(dterm/csum).lt.dcon).and.(abs(cosi(k,j+j+ixx))
     1             .gt.1.0e-01_knd))  go to 190
            dold=dnew
            j=j-1
180         continue
190       ce(k)=csum*dnuma*sgna
200       if((csum.eq.0.0e0_knd).and.fterm.ne.0.0e0_knd) naccc(k)=0
          if((csum*fterm).ne.0.0e0_knd) naccc(k)=ndec-2
     1                                    -log10(abs((fterm)/(csum)))
          if((csum.ne.0.0e0_knd).and.fterm.eq.0.0e0_knd) naccc(k)=ndec-2
          if((csum.eq.0.0e0_knd).and.fterm.eq.0.0e0_knd) naccc(k)=ndec
          if(naccc(k).gt.ndec-2) naccc(k)=ndec-2
          if(naccc(k).lt.0) naccc(k)=0
          if(naccc(k).gt.0) go to 210
          naccc(k)=0
          ce(k)=0.0e0_knd
          ced(k)=0.0e0_knd
          go to 280
c
c       compute the first derivative of the angular functions ce (when
c       iopang equals 2)
210       if(iopang.eq.1) go to 280
          if(abs(sine(k,3-ix)).gt.adec) go to 230
          ced(k)=0.0e0_knd
          jn=1
c          if(isq.eq.1.or.ix.eq.0) write(50,275) k,jfun,jn
c          if(isq.eq.-1.and.ix.eq.1) write(50,455) k,jfun,jn
          go to 280
230       doldd=1.0e0_knd
          fterm=0.0e0_knd
          cdsum=-real(l,knd)*sine(k,l+1)
          if(cdsum.gt.0.0e0_knd) fterm=cdsum
            do 240 j=jlow,lim
            coef=real(j+j+ix,knd)
            dnewd=doldd*enra(j)
            dterm=-dnewd*sine(k,j+j+ixx2)*coef
            if(dterm.gt.0.0e0_knd) fterm=fterm+dterm
            cdsum=cdsum+dterm
            if(cdsum.eq.0.0e0_knd) go to 235
            if((abs(dterm/cdsum).lt.dcon).and.(abs(sine(k,j+j+ixx2))
     1            .gt.1.0e-01_knd)) go to 250
235         doldd=dnewd
240         continue
250       if(j.gt.jangc) jangc=j
          jn=j
          doldd=1.0e0_knd
          j=l2
          ja=l2
          if(ja.eq.0) go to 270
            do 260 jj=1,ja
            coef=real(j+j+ix-2,knd)
            dnewd=doldd/enra(j)
            dterm=-dnewd*sine(k,j+j+ixx)*coef
            if(dterm.gt.0.0e0_knd) fterm=fterm+dterm
            cdsum=cdsum+dterm
            if((abs(dterm/cdsum).lt.dcon).and.(abs(sine(k,j+j+ixx))
     1            .gt.1.0e-01_knd)) go to 270
            doldd=dnewd
            j=j-1
260         continue
270       ced(k)=cdsum*dnuma*sgna
c          if(isq.eq.1.or.ix.eq.0) write(50,275) k,jfun,jn
c          if(isq.eq.-1.and.ix.eq.1) write(50,455) k,jfun,jn
c275       format(8x,'angle ',i3,'. ce converged in ',i6,' terms; ced ',
c     1           'converged in ',i6,' terms')
          if((cdsum.eq.0.0e0_knd).and.fterm.ne.0.0e0_knd) nacd=0
          if((cdsum*fterm).ne.0.0e0_knd) nacd=ndec-2
     1                                 -log10(abs((fterm)/(cdsum)))
          if((cdsum.ne.0.0e0_knd).and.fterm.eq.0.0e0_knd) nacd=ndec-2
          if((cdsum.eq.0.0e0_knd).and.fterm.eq.0.0e0_knd) nacd=ndec
          naccc(k)=min(naccc(k),nacd)
          if(naccc(k).lt.0) naccc(k)=0
          if(naccc(k).gt.0) go to 280
          ced(k)=0.0e0_knd
280       if(l.ne.0) go to 290
          se(k)=0.0e0_knd
          sed(k)=0.0e0_knd
          naccs(k)=ndec
290       continue
        jangs=0
        if(l.eq.0) go to 500
c
c  compute dnumb, the magnitude of the b coefficient with index = l
        dnumb=1.0e0_knd
        coef=1.0e0_knd
        jlow=l+2
        jterm=l2
          do 310 j=jlow,lims1,2
          jterm=jterm+1
          coef=coef*enrb(jterm)*enrb(jterm)
          dnumb=dnumb+coef
          if(abs(coef/dnumb).lt.dcon) go to 320
310       continue
320     jlow=l
        if(ix.eq.0) jlow=jlow-2
        jn=jterm
        if(l.lt.3) go to 340
        coef=1.0e0_knd
        jterm=l2
        j=jlow
          do 330 jj=2,jlow,2
          coef=coef/(enrb(jterm)*enrb(jterm))
          jterm=jterm-1
          j=j-2
          dnumb=dnumb+coef
          if(abs(coef/dnumb).lt.dcon) go to 340
330       continue
340     dnumb=1.0e0_knd/sqrt(dnumb)
        bsubl=dnumb*sgnb
c        if(isq.eq.1.or.ix.eq.0) write(50,345) jn,lim
c        if(isq.eq.-1.and.ix.eq.1) write(50,145) jn,lim
c345     format(5x,'se normalization series converged in ',
c     1         i6,' terms; ',i6,' terms available')
        jlow=l2+1
        jangs=0
           do 490 k=1,narg
c
c  compute the angular function se
          if(abs(sine(k,3-ix)).gt.adec) go to 350
          naccs(k)=ndec
          se(k)=0.0e0_knd
          jfun=1
c          if(iopang.eq.1.and.(isq.eq.1.or.ix.eq.0)) write(50,375) k,jfun
c          if(iopang.eq.1.and.isq.eq.-1.and.ix.eq.1) write(50,175) k,jfun
          go to 410
350       dold=1.0e0_knd
          ssum=sine(k,l+1)
          fterm=0.0e0_knd
          if(ssum.gt.0.0e0_knd) fterm=ssum
            do 360 j=jlow,lim
            dnew=dold*enrb(j)
            dterm=dnew*sine(k,j+j+ixx2)
            ssum=ssum+dterm
            if(dterm.gt.0.0e0_knd) fterm=fterm+dterm
            if((abs(dterm/ssum).lt.dcon).and.(abs(sine(k,j+j+ixx2))
     1            .gt.1.0e-01_knd)) go to 370
            dold=dnew
360         continue
370       if(j.gt.jangs) jangs=j
          jfun=j
c          if(iopang.eq.1.and.(isq.eq.1.or.ix.eq.0)) write(50,375) k,jfun
c          if(iopang.eq.1.and.isq.eq.-1.and.ix.eq.1) write(50,175) k,jfun
c375       format(8x,'angle ',i3,'. se converged in ',i6,' terms')
          jjupp=l2
          if(ix.eq.0) jjupp=l2-1
          if(jjupp.lt.1) go to 390
          dold=1.0e0_knd
          j=l2
            do 380 jj=1,jjupp
            dnew=dold/enrb(j)
            dterm=dnew*sine(k,j+j+ixx)
            ssum=ssum+dterm
            if(dterm.gt.0.0e0_knd) fterm=fterm+dterm
            if((abs(dterm/ssum).lt.dcon).and.(abs(sine(k,j+j+ixx))
     1            .gt.0.1e0_knd)) go to 390
            dold=dnew
            j=j-1
380         continue
390       se(k)=ssum*dnumb*sgnb
          if((ssum.eq.0.0e0_knd).and.fterm.ne.0.0e0_knd) naccs(k)=0
          if((ssum.ne.0.0e0_knd).and.fterm.ne.0.0e0_knd) naccs(k)=ndec-2
     1                                 -log10(abs((fterm)/(ssum)))
          if((ssum.eq.0.0e0_knd).and.fterm.eq.0.0e0_knd) naccs(k)=ndec
          if((ssum.ne.0.0e0_knd).and.fterm.eq.0.0e0_knd) naccs(k)=ndec-2
          if(naccs(k).gt.ndec-2) naccs(k)=ndec-2
          if(naccs(k).lt.0) naccs(k)=0
          if(naccs(k).gt.0) go to 410
          naccs(k)=0
          se(k)=0.0e0_knd
          sed(k)=0.0e0_knd
          go to 490
c
c       compute the first derivative of the angular functions se (when
c       iopang equals 2)
410       if(iopang.eq.1) go to 490
          if(abs(cosi(k,3-ix)).gt.adec) go to 415
          sed(k)=0.0e0_knd
          jn=1
c          if(isq.eq.1.or.ix.eq.0) write(50,455) k,jfun,jn
c          if(isq.eq.-1.and.ix.eq.1) write(50,275) k,jfun,jn
          go to 490
415       doldd=1.0e0_knd
          fterm=0.0e0_knd
          sdsum=real(l,knd)*(cosi(k,l+1))
          if(sdsum.gt.0.0e0_knd) fterm=sdsum
            do 440 j=jlow,lim
            coef=real(j+j+ix,knd)
            dnewd=doldd*enrb(j)
            dterm=dnewd*cosi(k,j+j+ixx2)*coef
            if(dterm.gt.0.0e0_knd) fterm=fterm+dterm
            sdsum=sdsum+dterm
            if(sdsum.eq.0.0e0_knd) go to 420
            if((abs(dterm/sdsum).lt.dcon).and.(abs(cosi(k,j+j+ixx2))
     1            .gt.1.0e-01_knd)) go to 450
420         doldd=dnewd
440         continue
450       if(j.gt.jangs) jangs=j
c          if(isq.eq.1.or.ix.eq.0) write(50,455) k,jfun,j
c          if(isq.eq.-1.and.ix.eq.1) write(50,275) k,jfun,j
c455       format(8x,'angle ',i3,'. se converged in ',i6,' terms; sed ',
c     1           'converged in ',i6,' terms')
          jjupp=l2
          if(ix.eq.0) jjupp=l2-1
          if(jjupp.lt.1) go to 470
          doldd=1.0e0_knd
          j=l2
            do 460 jj=1,jjupp
            coef=real(j+j+ix-2,knd)
            dnewd=doldd/enrb(j)
            dterm=dnewd*cosi(k,j+j+ixx)*coef
            if(dterm.gt.0.0e0_knd) fterm=fterm+dterm
            sdsum=sdsum+dterm
            if((abs(dterm/sdsum).lt.dcon).and.(abs(cosi(k,j+j+ixx))
     1            .gt.1.0e-01_knd)) go to 470
            doldd=dnewd
            j=j-1
460         continue
470       sed(k)=sdsum*dnumb*sgnb
          if((sdsum.eq.0.0e0_knd).and.fterm.ne.0.0e0_knd) nacd=0
          if((sdsum.ne.0.0e0_knd).and.fterm.ne.0.0e0_knd) nacd=ndec-2
     1                                 -log10(abs((fterm)/(sdsum)))
          if((sdsum.ne.0.0e0_knd).and.fterm.eq.0.0e0_knd) nacd=ndec-2
          if((sdsum.eq.0.0e0_knd).and.fterm.eq.0.0e0_knd) nacd=ndec
          naccs(k)=min(naccs(k),nacd)
          if(naccs(k).lt.0) naccs(k)=0
          if(naccs(k).gt.0) go to 490
          sed(k)=0.0e0_knd
490       continue
500       if(isq.eq.1) go to 600
            do 520 k=1,narg
            if(ix.ne.0) go to 510
            se(k)=-se(k)
            sed(k)=-sed(k)
            if(lx2.eq.0) go to 520
            ce(k)=-ce(k)
            ced(k)=-ced(k)
            se(k)=-se(k)
            sed(k)=-sed(k)
            go to 520
510         temp=ce(k)
            tempd=ced(k)
            ce(k)=se(k)
            ced(k)=sed(k)
            se(k)=temp
            sed(k)=tempd
            if(lx2.eq.0) go to 520
            ce(k)=-ce(k)
            ced(k)=-ced(k)
            se(k)=-se(k)
            sed(k)=-sed(k)
520         continue
600     return
        end
c
c
        subroutine m1bes (iopcs,l,cm,x1,isq,limr1,ndec,maxd,enra,
     1                    enrb,maxj,maxlp,nex,cbesf,cbesn,ibese,cbesdf,
     2                    cbesdr,m1bot,a01,ia01,b12,ib12,jbes,iflag,
     3                    nsubcs,m1c,im1e,m1dc,im1de)
c
c  purpose    : to calculate the Mathieu radial function of the
c               first kind and its first derivative with respect
c               to the traditional radial coordinate z, using an
c               expansion of cylindrical Bessel functions of the
c               first kind with argument c*sqrt(xi*xi-1)
c
c  parameters :
c
c     input   : iopcs  : =1 if cosine radial functions calculated
c                      : =2 if sine radial functions calculated
c               l      : l
c               cm     : magnitude of c
c               x1     : xi-1
c               isq    : integer = +1 if q is positive, = -1 if q is
c                        negative
c               limr1  : twice the maximum number of terms available
c               ndec   : number of decimal digits available in
c                        real(knd) arithmetic
c               maxd   : size of enr vector
c               enra   : a coefficient ratios
c               enrb   : b coefficient ratios
c               maxj   : size of cbesf vector
c               maxlp  : maximum  l value desired; size
c                        of the cbesn and ibese vectors
c               nex    : largest integer exponent available in
c                        real(knd) arithmetic
c               cbesf  : vector of ratios of cylindrical Bessel
c                        functions of the first kind
c               cbesn  : vector of characteristics of the Bessel
c                        functions
c               ibese  : vector of exponents of the Bessel functions
c               cbesdf : vector of ratios of first derivatives of the
c                        cylindrical Bessel functions
c               cbesdr : vector of ratios of first derivatives of the
c                        cylindrical Bessel functions to the
c                        corresponding Bessel functions
c               m1bot  : denominator sum for m1 and m1d that is
c                        calculated in dnorma or dnormb
c               a01    : characteristic of the first a coefficient
c               ia01   : exponent of the first a coefficient
c               a01    : characteristic of the first b coefficient
c               ia01   : exponent of the first b coefficient
c
c     output  : jbes   : maximum number of terms taken
c               iflag  : integer = 1 when backward series not used
c                        to calculate functions and first derivatives
c               nsubcs : larger of the subtraction errors encountered
c                        in the numerator series for m1 and m1d
c               m1c    : characteristic of Mathieu radial function
c                        of the first kind
c               im1e   : exponent of Mathieu radial function of the
c                        first kind [note that when q is
c                        negative (isq=-1) and l is odd, the Mathieu
c                        radial function m1 = i times m1c times
c                        10**(im1e), where i is the square root of -1.]
c               m1dc   : characteristic of derivative with respect
c                        to x of Mathieu radial function of the first
c                        kind
c               im1de  : exponent of derivative with respect to x of
c                        Mathieu radial function of the first kind
c                        [note that when qu is negative (isq=-1)
c                        and l is odd, the Mathieu radial function
c                        m1d = i times m1dc times 10**(im1de), where i
c                        is the square root of -1.]
c
        use param
c
c  real(knd) scalars and vectors
        real(knd) a1ol,a01,b12,cm,con1,dcon,dsqu,dnew,dnewd,dold,doldd,
     1            m1bot,m1c,m1dc,m1dpos,m1dstore,m1dtemp,m1dtempd,
     2            m1pos,m1temp,m1tempd,m1top,sqx2m1,term,termd,
     3            teste,testeo,xi,x1,x2m1
        real(knd) enr(maxd),enra(maxd),enrb(maxd),cbesdf(maxj),
     1            cbesdr(maxj),cbesf(maxj),cbesn(maxlp)
c
c  integer vector
        dimension ibese(maxlp)
c
c  convergence ratio dcon is set according to the number of decimal
c  digits available
        dcon=10.0e0_knd**(-ndec-2)
        xi=x1+1.0e0_knd
        x2m1=x1*(x1+2.0e0_knd)
        sqx2m1=sqrt(x2m1)
        l2=l/2
c  ix=0 for l even, ix=1 for l odd
        ix=l-2*l2
c  lx2=0 for l2 even, lx2=1 for l2 odd
        lx2=l2-2*(l2/2)
c  ix2=0 for (l+1)/2 even, ix2=1 for (l+1)/2 odd
        ix2=(l+1)/2-2*((l+1)/4)
        ics=ix
        if(iopcs.eq.2) ics=iabs(ix-1)
        dsqu=real(isq,knd)
        iflagd=0
        if(xi.lt.1.1e0_knd.and.iopcs.eq.1.and.ix.eq.1) iflagd=1
        if(x1.eq.0.0e0_knd) go to 220
        jlim=1
        if(ix.eq.0.and.iopcs.eq.2) jlim=2
        lim=limr1/2-ix
        con1=xi/sqx2m1
        nfac=nex-ndec
        teste=10.0e0_knd**nfac
        testeo=1.0e0_knd/teste
        im1dtemp=0
        im1temp=0
        if(isq.eq.-1) go to 25
        if(iopcs.eq.2) go to 15
          do 10 i=1,lim
          enr(i)=enra(i)
10        continue
        go to 50
15        do 20 i=1,lim
          enr(i)=enrb(i)
20        continue
        go to 50
25      if(ics.eq.1) go to 35
          do 30 i=1,lim
          enr(i)=enra(i)
30        continue
        go to 50
35        do 40 i=1,lim
          enr(i)=enrb(i)
40        continue
50      continue
c
c  compute ratio of 1st coeffficient to lth coefficient when l
c  is odd, xi is less than 1.1 but not equal to 1, and cosine functions
c  are being computed
c
        if(iflagd.eq.0) go to 70
          a1ol=1.0e0_knd
          ia1ol=0
          if(l.eq.1) go to 70
          do 60 i=1,l2
          a1ol=a1ol/enr(i)
          iterm=int(log10(abs(a1ol)))
          a1ol=a1ol*(10.0e0_knd**(-iterm))
          ia1ol=ia1ol+iterm
60        continue
70      continue
c
c  compute radial function of the first kind m1
c
c  forward summation of numerator series for m1 and m1d
        m1tempd=0.0e0_knd
        m1dtempd=0.0e0_knd
        dold=1.0e0_knd
        doldd=1.0e0_knd
        m1temp=dold
        m1dtemp=doldd
        jlow=l2+1
        if(ics.eq.0) go to 80
        m1temp=real(l,knd)*dold
        m1dtemp=real(l,knd)*doldd
80      m1pos=m1temp
        m1dpos=m1dtemp
        jbes=l2
        if(iflag.eq.1) go to 120
          do 100 j=jlow,lim
          jj=j+j+ix
          dnew=dold*enr(j)*cbesf(jj)
          dnewd=doldd*enr(j)*cbesdf(jj)
          term=dnew
          termd=dnewd
          if(ics.eq.0) go to 90
          term=dnew*real(jj,knd)
          termd=dnewd*real(jj,knd)
90        if(term.gt.0.0e0_knd) m1pos=m1pos+term
          if(termd.gt.0.0e0_knd) m1dpos=m1dpos+termd
          m1temp=m1temp+term
          m1dtemp=m1dtemp+termd
          if((abs(term/m1temp)+abs(termd/m1dtemp)).lt.dcon) go to 110
          dold=dnew
          doldd=dnewd
100       continue
110     jbes=j
        if(iflagd.eq.0.or.l.ne.1) go to 120
        m1tempd=m1temp-1.0e0_knd
        m1dtempd=m1dtemp-1.0e0_knd
120	continue
c
c  backward summation of numerator series for m1
        if (l2.lt.jlim) go to 175
        dold=1.0e0_knd
        doldd=1.0e0_knd
          do 160 j=l2,jlim,-1
          jj=j+j+ix
          dnew=dold/(enr(j)*cbesf(jj))
          dnewd=doldd/(enr(j)*cbesdf(jj))
          term=dnew
          termd=dnewd
          if(ics.eq.0) go to 140
          term=dnew*real(jj-2,knd)
          termd=dnewd*real(jj-2,knd)
140       if(term.gt.0.0e0_knd) m1pos=m1pos+term
          if(termd.gt.0.0e0_knd) m1dpos=m1dpos+termd
          m1temp=m1temp+term
          m1dtemp=m1dtemp+termd
          if((abs(term/m1temp)+abs(termd/m1dtemp)).lt.dcon) go to 170
            if(abs(m1temp).gt.teste.or.abs(m1dtemp).gt.teste) then
            m1temp=m1temp*testeo
            m1dtemp=m1dtemp*testeo
            dnew=dnew*testeo
            dnewd=dnewd*testeo
            m1pos=max(m1pos*testeo,dcon*abs(m1temp))
            m1dpos=max(m1dpos*testeo,dcon*abs(m1dtemp))
            im1temp=im1temp+nfac
            end if
150       dold=dnew
          doldd=dnewd
160       continue
170     if(jj.ne.3) iflagd=0
        im1dtemp=im1temp
        if(iflagd.eq.0) go to 175
        m1tempd=m1temp-dnew
        m1dtempd=m1dtemp-dnewd
175     m1top=m1temp*cbesn(l+1)
        im1tope=im1temp+ibese(l+1)
        iterm=0
        if(m1top.ne.0.0e0_knd) iterm=log10(abs(m1top))
        im1tope=im1tope+iterm
        m1top=m1top*10.0e0_knd**(-iterm)
        nsub=0
        if((m1temp*m1pos).ne.0.0e0_knd) nsub=
     1                             int(log10(abs(m1pos/m1temp)))
        if(nsub.lt.0) nsub=0
        if(nsub.gt.ndec) nsub=ndec
        if(m1temp.eq.0.0e0_knd.and.m1pos.ne.0.0e0_knd) nsub=ndec
        nsubd=0
        if((m1dtemp*m1dpos).ne.0.0e0_knd) nsubd=
     1                             int(log10(abs(m1dpos/m1dtemp)))
        if(nsubd.lt.0) nsubd=0
        if(nsubd.gt.ndec) nsubd=ndec
        if(m1dtemp.eq.0.0e0_knd.and.m1dpos.ne.0.0e0_knd) nsubd=ndec
        nsubcs=max(nsub,nsubd)
c        if(iopcs.eq.1.and.iflag.eq.0) write(40,180) jbes,lim,nsubcs
c180     format(3x,'mc1/mc1d numerator converged in ',i5,' terms; ',i5,
c     1         ' available; sub. error = ',i3,' digits.')
c        if(iopcs.eq.2.and.iflag.eq.0) write(40,190) jbes,lim,nsubcs
c190     format(3x,'ms1/ms1d numerator converged in ',i5,' terms; ',i5,
c     1         ' available; sub. error = ',i3,' digits.')
c        if(iopcs.eq.1.and.iflag.eq.1) write(40,200) jbes,lim,nsubcs
c200     format(3x,'mc1/mc1d numerator converged in ',i5,' terms; ',i5,
c     1         ' available; fwd series not used; sub. error = ',i3,
c     2         ' digits.')
c        if(iopcs.eq.2.and.iflag.eq.1) write(40,210) jbes,lim,nsubcs
c210     format(3x,'ms1/ms1d numerator converged in ',i5,' terms; ',i5,
c     1         ' available; fwd series not used; sub. error = ',i3,
c     2         ' digits.')
c
c  combining numerator with denominator series computed in dnorma for
c  mc1 or in dnormb for ms1
        m1c=m1top/m1bot
        if(ics.eq.1) m1c=m1c*con1
        go to 240
220     if(isq.eq.-1) go to 230
        m1c=0.0e0_knd
        im1e=0
        if(iopcs.eq.1) m1c=a01/m1bot
        if(iopcs.eq.1) im1e=ia01
        if(ix.eq.1) m1c=0.5e0_knd*m1c*cm
        go to 240
230     m1c=0.0e0_knd
        im1e=0
        if(iopcs.eq.1.and.ix.eq.0) m1c=a01/m1bot
        if(iopcs.eq.1.and.ix.eq.0) im1e=ia01
        if(iopcs.eq.1.and.ix.eq.1) m1c=0.5e0_knd*cm*b12/m1bot
        if(iopcs.eq.1.and.ix.eq.1) im1e=ib12
240     if(isq.eq.1.and.iopcs.eq.1.and.ix2.eq.1) m1c=-m1c
        if(isq.eq.1.and.iopcs.eq.2.and.lx2.eq.1) m1c=-m1c
        if(isq.eq.-1.and.iopcs.eq.1.and.lx2.eq.1) m1c=-m1c
        if(isq.eq.-1.and.iopcs.eq.2.and.lx2.eq.1) m1c=-m1c
        iterm=0
        if(m1c.ne.0.0e0_knd) iterm=int(log10(abs(m1c)))
        if(x1.eq.0.0e0_knd) im1e=im1e+iterm
        if(x1.ne.0.0e0_knd) im1e=im1tope+iterm
        m1c=m1c*10.0e0_knd**(-iterm)
        if(abs(m1c).ge.1.0e0_knd.or.m1c.eq.0.0e0_knd) go to 250
        m1c=m1c*10.0e0_knd
        im1e=im1e-1
250	continue
        if(x1.eq.0.0e0_knd) go to 270
        if(iflagd.eq.1) m1temp=m1tempd
        if(iflagd.eq.1) m1dtemp=m1dtempd
        m1dtemp=cm*xi*m1dtemp*cbesdr(l+1)
        m1dstore=m1dtemp
        if(ics.eq.1) m1dtemp=m1dtemp-m1temp/(xi*sqx2m1)
        idsub=-log10(abs(m1dtemp/m1dstore))
        if(idsub.lt.0) idsub=0
        m1dc=m1dtemp*cbesn(l+1)/m1bot
        if(iflagd.eq.0) go to 255
        term=x2m1*cbesdr(2)*cbesn(2)
        termd=term-dsqu*cbesn(3)*(10.0e0_knd**(ibese(3)-ibese(2)))
        idsub1=-log10(abs(termd/term))
        if(idsub1.lt.0) idsub1=0
        termd=termd*a1ol*(cm/xi)*(10.0e0_knd**(ia1ol+ibese(2)-
     1        ibese(l+1)-im1dtemp))/m1bot
        m1dc=m1dc+termd
        idsub2=-log10(abs(m1dc/termd))
        if(idsub2.lt.0) idsub2=0
        idsub=idsub+idsub1+idsub2
255     continue
c        if(idsub.ne.0) write(40,260) idsub
c260     format(24x,'subtraction error in forming m1d =',i3,' digits.')
        if(ics.eq.1) m1dc=m1dc*con1
        go to 290
270     if(isq.eq.-1) go to 280
        m1dc=0.0e0_knd
        im1de=0
        if(iopcs.eq.2) m1dc=0.5e0_knd*cm*b12/m1bot
        if(iopcs.eq.2) im1de=ib12
        if(ix.eq.0) m1dc=0.5e0_knd*cm*m1dc
        go to 290
280     m1dc=0.0e0_knd
        im1de=0
        if(iopcs.eq.2.and.ix.eq.0) m1dc=0.25e0_knd*cm*cm*b12/m1bot
        if(iopcs.eq.2.and.ix.eq.0) im1de=ib12
        if(iopcs.eq.2.and.ix.eq.1) m1dc=0.5e0_knd*cm*a01/m1bot
        if(iopcs.eq.2.and.ix.eq.1) im1de=ia01
290     if(isq.eq.1.and.iopcs.eq.1.and.ix2.eq.1) m1dc=-m1dc
        if(isq.eq.1.and.iopcs.eq.2.and.lx2.eq.1) m1dc=-m1dc
        if(isq.eq.-1.and.iopcs.eq.1.and.lx2.eq.1) m1dc=-m1dc
        if(isq.eq.-1.and.iopcs.eq.2.and.lx2.eq.1) m1dc=-m1dc
        iterm=0
        if(m1dc.ne.0.0e0_knd) iterm=int(log10(abs(m1dc)))
        if(x1.ne.0.0e0_knd) im1de=iterm+ibese(l+1)+im1dtemp
        if(x1.eq.0.0e0_knd) im1de=im1de+iterm
        m1dc=m1dc*10.0e0_knd**(-iterm)
        if(abs(m1dc).ge.1.0e0_knd.or.m1dc.eq.0.0e0_knd) go to 300
        m1dc=10.0e0_knd*m1dc
        im1de=im1de-1
300     mfac=im1tope-ibese(l+1)
        if(mfac.gt.(ndec+5)) iflag=1
        if(mfac.le.(ndec+5)) iflag=0
        return
        end
c
c
        subroutine m1bpe (iopcs,iss,ismax,l,cm,x1,limbpe,ndec,limd,maxd,
     1                    enra,enrb,maxj,maxlp,nex,cbesf1,cbesn1,ibese1,
     2                    cbesdf1,cbesdr1,cbesf,cbesn,ibese,cbesdf,
     3                    cbesdr,a01,ia01,b12,ib12,jbpe,nsub,m1c,im1e,
     4                    m1dc,im1de)
c
c  purpose    : to calculate the Mathieu radial function of the
c               first kind and its first derivative with respect
c               to the traditional radial function z when q is
c               negative, using a series of products of modified
c               cylindrical Bessel functions of the first kind,
c               with the integer order offset s (is) of the Bessel
c               function for input as a free parameter
c
c  parameters :
c
c     input   : iopcs  : =1 if cosine radial functions calculated
c                      : =2 if sine radial functions calculated
c               iss    : integer offset minus l from calculation
c                        for radial function of the same kind
c                        (sine or cosine) for previous l
c               ismax  : maximum allowable value for the integer offset
c               l      : l
c               cm     : magnitude of c (c being positive imaginary)
c               x1     : xi-1
c               limbpe : twice the maximum number of terms available
c               ndec   : number of decimal digits available in
c                        real(knd) arithmetic
c               limd   : number of a and b coefficient ratios that have
c                        been calculated for this case
c               maxd   : dimension of enra and enrb vectors
c               enra   : vector of a coefficient ratios
c               enrb   : vector of b coefficient ratios
c               maxj   : dimension of cbesf, cbesdf, cbesdr, cbesf1,
c                        cbesdf1 and cbesdr1 vectors
c               maxlp  : dimension of cbesn, ibese, cbesn1, and
c                        ibese1 vectors
c               nex    : largest integer exponent abvailable in
c                        real(knd) arithmetic
c               cbesf1 : vector of ratios of cylindrical Bessel
c                        functions of the first kind with argument
c                        c*[xi - sqrt(xi*xi - 1)]/2
c               cbesn1 : vector of characteristics of these Bessel
c                        functions
c               ibese1 : vector of exponents of these Bessel functions
c               cbesdf1: vector of ratios of first derivatives of
c                        these Bessel functions
c               cbesdr1: vector of ratios of first derivatives of
c                        these Bessel functions to the corresponding
c                        Bessel functions
c               cbesf  : vector of ratios of cylindrical Bessel
c                        functions of the first kind with argument
c                        c*[xi + sqrt(xi*xi - 1)]/2
c               cbesn  : vector of characteristics of these Bessel
c                        functions
c               ibese  : vector of exponents of these Bessel functions
c               cbesdf : vector of ratios of first derivatives of these
c                        Bessel functions
c               cbesdr : vector of ratios of first derivatives of
c                        these Bessel functions to the corresponding
c                        Bessel functions
c               a01    : characteristic of the first a coefficient
c               ia01   : exponent of the first a coefficient
c               b12    : characteristic of the first b coefficient
c               ib12   : exponent of the first b coefficient
c
c     output  : jbpe   : twice the maximum number of terms taken
c               nsub   : maximum number of digits of subtraction
c                        error encountered in the calculation of
c                        the function and its first derivative
c               m1c    : characteristic of Mathieu radial function
c                        of the first kind.
c               im1e   : exponent of Mathieu radial function of the
c                        first kind
c               m1dc   : characteristic of derivative with respect
c                        to x of Mathieu radial function of the first
c                        kind
c               im1de  : exponent of derivative with respect to x of
c                        Mathieu radial function of the first kind
c
        use param
c
c  real(knd) scalars and vectors
        real(knd) a01,b12,cm,coef,coefa,coefb,dcon,den,dnewa,dnewb,
     1            dolda,doldas,doldb,doldbs,m1bot,m1c,m1dc,m1dpos,
     2            m1dtemp,m1pos,m1temp,sqx2m1,term,terma,termb,termd,
     3            termd1,termd1a,termd1b,termd2,termd2a,termd2b,
     4            teste,testeo,u1,u2,xi,x1,x2m1
        real(knd) enr(maxd),enra(maxd),enrb(maxd),cbesdf1(maxj),
     1            cbesdr1(maxj),cbesf1(maxj),cbesn1(maxlp),cbesdf(maxj),
     2            cbesdr(maxj),cbesf(maxj),cbesn(maxlp)
c
c  integer vectors
        dimension ibese1(maxlp),ibese(maxlp)
c
c  convergence ratio dcon is set according to the number of decimal
c  digits available
        dcon=10.0e0_knd**(-ndec-2)
        xi=x1+1.0e0_knd
        x2m1=x1*(x1+2.0e0_knd)
        sqx2m1=sqrt(x2m1)
        u1=0.5e0_knd*cm/(xi+sqx2m1)
        u2=0.5e0_knd*cm*(xi+sqx2m1)
        l2=l/2
c
c  ix=0 for l even, ix=1 for l odd
        ix=l-2*l2
        ics=ix
        if(iopcs.eq.2) ics=iabs(ix-1)
        if(iopcs.eq.1.and.ix.eq.0) iopt=10
        if(iopcs.eq.1.and.ix.eq.1) iopt=11
        if(iopcs.eq.2.and.ix.eq.0) iopt=20
        if(iopcs.eq.2.and.ix.eq.1) iopt=21
        lim=limbpe/2-ix
        nfac=nex-ndec
        teste=10.0e0_knd**nfac
        testeo=1.0e0_knd/teste
        im1dtemp=0
        im1temp=0
        joff=0
        if(iopt.eq.20) joff=1
        jlim=1+joff
        is=l+iss
        if(is.gt.ismax) is=ismax
        lims=is+joff+100+int(cm/10.0)
        limmax=max(lims,lim)
        if(limmax.gt.limd-1) limmax=limd-1
        if(ics.eq.1) go to 20
          do 10 i=1,limmax
          enr(i)=enra(i)
10        continue
          m1bot=a01
          ibot=ia01
          go to 40
20        do 30 i=1,limmax
          enr(i)=enrb(i)
30        continue
          m1bot=b12
          ibot=ib12
40      continue
        den=m1bot
        iden=ibot
          if(is.eq.0) go to 60
          do 50 k=1,is
          den=den*enr(k+joff)
          iterm=int(log10(abs(den)))
          den=den*(10.0e0_knd**(-iterm))
          iden=iden+iterm
50        continue
60      continue
c
c  compute radial function of the first kind m1
c
c    forward summation of product series for m1 and m1d
        coefa=cbesn1(-l2+is+1+joff)*cbesn(l2+is+1+joff)
        inda=ibese1(-l2+is+1+joff)+ibese(l2+is+1+joff)
        coefb=cbesn1(l2+is+1+joff)*cbesn(-l2+is+1+joff)
        indb=ibese1(l2+is+1+joff)+ibese(-l2+is+1+joff)
        if(inda.lt.indb) go to 70
        dolda=1.0e0_knd
        doldb=(coefb/coefa)*(10.0e0_knd**(indb-inda))
        coef=coefa
        icoef=inda
        go to 80
70      doldb=1.0e0_knd
        dolda=(coefa/coefb)*(10.0e0_knd**(inda-indb))
        coef=coefb
        icoef=indb
80      doldas=dolda
        doldbs=doldb
        if(iopt.eq.10) terma=1.0e0_knd
        if(iopt.eq.10) termb=1.0e0_knd
        if(iopt.eq.11) terma=cbesf(l2+is+1)
        if(iopt.eq.11) termb=cbesf1(l2+is+1)
        if(iopt.eq.21) terma=cbesf(l2+is+1)
        if(iopt.eq.21) termb=-cbesf1(l2+is+1)
        if(iopt.eq.20) terma=1.0e0_knd
        if(iopt.eq.20) termb=-1.0e0_knd
        m1temp=dolda*terma+doldb*termb
        if(iopt.eq.10) termd1a=cbesdr1(is-l2+1)
        if(iopt.eq.10) termd2a=cbesdr(l2+is+1)
        if(iopt.eq.10) termd1b=cbesdr1(l2+is+1)
        if(iopt.eq.10) termd2b=cbesdr(is-l2+1)
        if(iopt.eq.20) termd1a=cbesdr1(is-l2+2)
        if(iopt.eq.20) termd2a=cbesdr(l2+is+2)
        if(iopt.eq.20) termd1b=-cbesdr1(l2+is+2)
        if(iopt.eq.20) termd2b=-cbesdr(is-l2+2)
        if(iopt.eq.11) termd1a=cbesdr1(is-l2+1)*cbesf(l2+is+1)
        if(iopt.eq.11) termd2a=cbesdr(l2+is+2)*cbesf(l2+is+1)
        if(iopt.eq.11) termd1b=cbesdr1(l2+is+2)*cbesf1(l2+is+1)
        if(iopt.eq.11) termd2b=cbesdr(is-l2+1)*cbesf1(l2+is+1)
        if(iopt.eq.21) termd1a=cbesdr1(is-l2+1)*cbesf(l2+is+1)
        if(iopt.eq.21) termd2a=cbesdr(l2+is+2)*cbesf(l2+is+1)
        if(iopt.eq.21) termd1b=-cbesdr1(l2+is+2)*cbesf1(l2+is+1)
        if(iopt.eq.21) termd2b=-cbesdr(is-l2+1)*cbesf1(l2+is+1)
        termd1=dolda*termd1a+doldb*termd1b
        termd2=dolda*termd2a+doldb*termd2b
        m1dtemp=-u1*termd1+u2*termd2
        m1pos=0.0e0_knd
        if(m1temp.gt.0.0e0_knd) m1pos=m1temp
        m1dpos=0.0e0_knd
        if(m1dtemp.gt.0.0e0_knd) m1dpos=m1dtemp
        jbpe=l2
        lflag=0
c
c      forward summation for is => j
          do 100 j=l2+1,is+joff
          dnewa=-dolda*enr(j)*cbesf(j+is+joff)/cbesf1(is+joff-j+1)
          dnewb=-doldb*enr(j)*cbesf1(j+is+joff)/cbesf(is+joff-j+1)
          if(iopt.eq.10) terma=dnewa
          if(iopt.eq.10) termb=dnewb
          if(iopt.eq.11) terma=dnewa*cbesf(j+is+1)
          if(iopt.eq.11) termb=dnewb*cbesf1(j+is+1)
          if(iopt.eq.21) terma=dnewa*cbesf(j+is+1)
          if(iopt.eq.21) termb=-dnewb*cbesf1(j+is+1)
          if(iopt.eq.20) terma=dnewa
          if(iopt.eq.20) termb=-dnewb
          term=terma+termb
          if(iopt.eq.10) termd1a=dnewa*cbesdr1(is-j+1)
          if(iopt.eq.10) termd2a=dnewa*cbesdr(j+is+1)
          if(iopt.eq.10) termd1b=dnewb*cbesdr1(j+is+1)
          if(iopt.eq.10) termd2b=dnewb*cbesdr(is-j+1)
          if(iopt.eq.20) termd1a=dnewa*cbesdr1(is-j+2)
          if(iopt.eq.20) termd2a=dnewa*cbesdr(j+is+2)
          if(iopt.eq.20) termd1b=-dnewb*cbesdr1(j+is+2)
          if(iopt.eq.20) termd2b=-dnewb*cbesdr(is-j+2)
          if(iopt.eq.11) termd1a=dnewa*cbesdr1(is-j+1)*cbesf(j+is+1)
          if(iopt.eq.11) termd2a=dnewa*cbesdr(j+is+2)*cbesf(j+is+1)
          if(iopt.eq.11) termd1b=dnewb*cbesdr1(j+is+2)*cbesf1(j+is+1)
          if(iopt.eq.11) termd2b=dnewb*cbesdr(is-j+1)*cbesf1(j+is+1)
          if(iopt.eq.21) termd1a=dnewa*cbesdr1(is-j+1)*cbesf(j+is+1)
          if(iopt.eq.21) termd2a=dnewa*cbesdr(j+is+2)*cbesf(j+is+1)
          if(iopt.eq.21) termd1b=-dnewb*cbesdr1(j+is+2)*cbesf1(j+is+1)
          if(iopt.eq.21) termd2b=-dnewb*cbesdr(is-j+1)*cbesf1(j+is+1)
          termd1=termd1a+termd1b
          termd2=termd2a+termd2b
          termd=-u1*termd1+u2*termd2
          if(term.gt.0.0e0_knd) m1pos=m1pos+term
          if(termd.gt.0.0e0_knd) m1dpos=m1dpos+termd
          m1temp=m1temp+term
          m1dtemp=m1dtemp+termd
          if(m1temp*m1dtemp.eq.0.0e0_knd) go to 120
          if((abs(term/m1temp)+abs(termd/m1dtemp)).lt.dcon) go to 120
            if(abs(m1temp).gt.teste.or.abs(m1dtemp).gt.teste) then
            m1temp=m1temp*testeo
            m1dtemp=m1dtemp*testeo
            dnewa=dnewa*testeo
            dnewb=dnewb*testeo
            m1pos=max(m1pos*testeo,dcon*abs(m1temp))
            m1dpos=max(m1dpos*testeo,dcon*abs(m1dtemp))
            m1dpos=m1dpos*testeo
            im1temp=im1temp+nfac
            im1dtemp=im1dtemp+nfac
            lflag=1
            end if
90        dolda=dnewa
          doldb=dnewb
100       continue
c
c      forward summation for j > is
          do 110 j=is+joff+1,lim
          dnewa=-dolda*enr(j)*cbesf1(j-is-joff)*cbesf(j+is+joff)
          dnewb=-doldb*enr(j)*cbesf1(j+is+joff)*cbesf(j-is-joff)
          if(iopt.eq.10) terma=dnewa
          if(iopt.eq.10) termb=dnewb
          if(iopt.eq.11) terma=dnewa*cbesf(j+is+1)
          if(iopt.eq.11) termb=dnewb*cbesf1(j+is+1)
          if(iopt.eq.21) terma=dnewa*cbesf(j+is+1)
          if(iopt.eq.21) termb=-dnewb*cbesf1(j+is+1)
          if(iopt.eq.20) terma=dnewa
          if(iopt.eq.20) termb=-dnewb
          term=terma+termb
          if(iopt.eq.10) termd1a=dnewa*cbesdr1(j-is+1)
          if(iopt.eq.10) termd2a=dnewa*cbesdr(j+is+1)
          if(iopt.eq.10) termd1b=dnewb*cbesdr1(j+is+1)
          if(iopt.eq.10) termd2b=dnewb*cbesdr(j-is+1)
          if(iopt.eq.20) termd1a=dnewa*cbesdr1(j-is)
          if(iopt.eq.20) termd2a=dnewa*cbesdr(j+is+2)
          if(iopt.eq.20) termd1b=-dnewb*cbesdr1(j+is+2)
          if(iopt.eq.20) termd2b=-dnewb*cbesdr(j-is)
          if(iopt.eq.11) termd1a=dnewa*cbesdr1(j-is+1)*cbesf(j+is+1)
          if(iopt.eq.11) termd2a=dnewa*cbesdr(j+is+2)*cbesf(j+is+1)
          if(iopt.eq.11) termd1b=dnewb*cbesdr1(j+is+2)*cbesf1(j+is+1)
          if(iopt.eq.11) termd2b=dnewb*cbesdr(j-is+1)*cbesf1(j+is+1)
          if(iopt.eq.21) termd1a=dnewa*cbesdr1(j-is+1)*cbesf(j+is+1)
          if(iopt.eq.21) termd2a=dnewa*cbesdr(j+is+2)*cbesf(j+is+1)
          if(iopt.eq.21) termd1b=-dnewb*cbesdr1(j+is+2)*cbesf1(j+is+1)
          if(iopt.eq.21) termd2b=-dnewb*cbesdr(j-is+1)*cbesf1(j+is+1)
          termd1=termd1a+termd1b
          termd2=termd2a+termd2b
          termd=-u1*termd1+u2*termd2
          if(term.gt.0.0e0_knd) m1pos=m1pos+term
          if(termd.gt.0.0e0_knd) m1dpos=m1dpos+termd
          m1temp=m1temp+term
          m1dtemp=m1dtemp+termd
          if(m1temp*m1dtemp.eq.0.0e0_knd) go to 120
          if((abs(term/m1temp)+abs(termd/m1dtemp)).lt.dcon) go to 120
          dolda=dnewa
          doldb=dnewb
110       continue
120     jbpe=j
c
c  backward summation of product series for m1 and m1d
        if(l2.lt.jlim.or.lflag.eq.1) go to 140
        dolda=doldas
        doldb=doldbs
          do 130 j=l2,jlim,-1
          dnewa=-dolda*cbesf1(is-j+1+joff)/(enr(j)*cbesf(j+is+joff))
          dnewb=-doldb*cbesf(is-j+1+joff)/(enr(j)*cbesf1(j+is+joff))
          if(iopt.eq.10) terma=dnewa
          if(iopt.eq.10) termb=dnewb
          if(iopt.eq.11) terma=dnewa*cbesf(j+is)
          if(iopt.eq.11) termb=dnewb*cbesf1(j+is)
          if(iopt.eq.21) terma=dnewa*cbesf(j+is)
          if(iopt.eq.21) termb=-dnewb*cbesf1(j+is)
          if(iopt.eq.20) terma=dnewa
          if(iopt.eq.20) termb=-dnewb
          term=terma+termb
          if(iopt.eq.10) termd1a=dnewa*cbesdr1(is-j+2)
          if(iopt.eq.10) termd2a=dnewa*cbesdr(j+is)
          if(iopt.eq.10) termd1b=dnewb*cbesdr1(j+is)
          if(iopt.eq.10) termd2b=dnewb*cbesdr(is-j+2)
          if(iopt.eq.20) termd1a=dnewa*cbesdr1(is-j+3)
          if(iopt.eq.20) termd2a=dnewa*cbesdr(j+is+1)
          if(iopt.eq.20) termd1b=-dnewb*cbesdr1(j+is+1)
          if(iopt.eq.20) termd2b=-dnewb*cbesdr(is-j+3)
          if(iopt.eq.11) termd1a=dnewa*cbesdr1(is-j+2)*cbesf(j+is)
          if(iopt.eq.11) termd2a=dnewa*cbesdr(j+is+1)*cbesf(j+is)
          if(iopt.eq.11) termd1b=dnewb*cbesdr1(j+is+1)*cbesf1(j+is)
          if(iopt.eq.11) termd2b=dnewb*cbesdr(is-j+2)*cbesf1(j+is)
          if(iopt.eq.21) termd1a=dnewa*cbesdr1(is-j+2)*cbesf(j+is)
          if(iopt.eq.21) termd2a=dnewa*cbesdr(j+is+1)*cbesf(j+is)
          if(iopt.eq.21) termd1b=-dnewb*cbesdr1(j+is+1)*cbesf1(j+is)
          if(iopt.eq.21) termd2b=-dnewb*cbesdr(is-j+2)*cbesf1(j+is)
          termd1=termd1a+termd1b
          termd2=termd2a+termd2b
          termd=-u1*termd1+u2*termd2
          if(term.gt.0.0e0_knd) m1pos=m1pos+term
          if(termd.gt.0.0e0_knd) m1dpos=m1dpos+termd
          m1temp=m1temp+term
          m1dtemp=m1dtemp+termd
          if(m1temp*m1dtemp.eq.0.0e0_knd) go to 140
          if((abs(term/m1temp)+abs(termd/m1dtemp)).lt.dcon) go to 140
          dolda=dnewa
          doldb=dnewb
130       continue
140     nsub=0
        if((m1temp*m1pos).ne.0.0e0_knd) nsub=
     1                             int(log10(abs(m1pos/m1temp)))
        if(nsub.lt.0) nsub=0
        if(nsub.gt.ndec) nsub=ndec
        if(m1temp.eq.0.0e0_knd.and.m1pos.ne.0.0e0_knd) nsub=ndec
        nsubd=0
        if((m1dtemp*m1dpos).ne.0.0e0_knd) nsubd=
     1                             int(log10(abs(m1dpos/m1dtemp)))
        if(nsubd.lt.0) nsubd=0
        if(nsubd.gt.ndec) nsubd=ndec
        if(m1dtemp.eq.0.0e0_knd.and.m1dpos.ne.0.0e0_knd) nsubd=ndec
c        if(iopcs.eq.1) write(40,150) jbpe,lim,nsub,nsubd,is
c150     format(3x,'mc1/mc1d converged in ',i5,' of ',i5,' terms;',
c     1         ' sub. errors =',i3,' and',i3,' digits. s = ',i6)
c        if(iopcs.eq.2) write(40,160) jbpe,lim,nsub,nsubd,is
c160     format(3x,'ms1/ms1d converged in ',i5,' of ',i5,' terms;',
c     1         ' sub. errors =',i3,' and',i3,' digits. s = ',i6)
          if(nsub.ne.0.or.nsubd.ne.0) then
          issinc=(nsub+nsubd)*(1.5-min(0,int(3.0e0_knd*(log10(x1)))))+
     1            max(0,2*int(log10(cm)))
          if(l.lt.200.and.cm.gt.1000.0e0_knd.and.x1.lt.0.1e0_knd)
     1          issinc=issinc*cm/1000.0e0_knd
          iss=iss+issinc
          end if
        iterm=log10(abs(m1temp))
        m1temp=m1temp*10.0e0_knd**(-iterm)
        m1c=m1temp*coef/den
        if(is.eq.0) m1c=m1c/2.0e0_knd
        im1e=im1temp+iterm+icoef-iden
        iterm=log10(abs(m1c))
        im1e=im1e+iterm
        m1c=m1c*10.0e0_knd**(-iterm)
        if(abs(m1c).ge.1.0e0_knd) go to 170
        m1c=m1c*10.0e0_knd
        im1e=im1e-1
170	continue
        m1c=-m1c
        if(iopt.eq.10) m1c=-m1c
        if(ix.eq.1) m1c=-m1c
        if(2*(is/2).ne.is) m1c=-m1c
        iterm=log10(abs(m1dtemp))
        m1dtemp=m1dtemp*10.0e0_knd**(-iterm)
        m1dc=m1dtemp*coef/den
        if(is.eq.0) m1dc=m1dc/2.0e0_knd
        im1de=im1dtemp+iterm+icoef-iden
        iterm=log10(abs(m1dc))
        im1de=im1de+iterm
        m1dc=m1dc*10.0e0_knd**(-iterm)
        if(abs(m1dc).ge.1.0e0_knd) go to 180
        m1dc=m1dc*10.0e0_knd
        im1de=im1de-1
180     continue
        m1dc=-m1dc
        if(iopt.eq.10) m1dc=-m1dc
        if(ix.eq.1) m1dc=-m1dc
        if(2*(is/2).ne.is) m1dc=-m1dc
        return
        end
c
c
        subroutine m1bpe0 (iopcs,iss,ismax,l,cm,limbpe,ndec,limd,
     1                     maxd,enra,enrb,maxj,maxlp,nex,cbesf,cbesn,
     2                     ibese,cbesdf,cbesdr,a01,ia01,b12,ib12,jbpe,
     3                     nsub,m1c,im1e,m1dc,im1de)
c
c  purpose    : to calculate the Mathieu radial function of the
c               first kind and its first derivative with respect
c               to the traditional radial function z when q is negative
c               and xi = 0, using a series of products of modified
c               cylindrical Bessel functions of the first kind, with
c               the integer order offset s (is) of the Bessel function
c               for input as a free parameter
c
c  parameters :
c
c     input   : iopcs  : =1 if cosine radial functions calculated
c                      : =2 if sine radial functions calculated
c               iss    : integer offset minus l from calculation
c                        for radial function of the same kind
c                        (sine or cosine) for previous l
c               ismax  : maxable allowable value for the integer offset
c               l      : l
c               cm     : magnitude of c (c being positive imaginary)
c               limbpe : twice the maximum number of terms available
c               ndec   : number of decimal digits available in
c                        real(knd) arithmetic
c               limd   : number of a and b coefficient ratios that have
c                        been calculated for this case
c               maxd   : dimension of enra and enrb vectors
c               enra   : vector of a coefficient ratios
c               enrb   : vector of b coefficient ratios
c               maxj   : dimension of cbesf, cbesdf, and cbesdr vectors
c               maxlp  : dimension of cbesn and ibese vectors
c               nex    : largest integer exponent abvailable in
c                        real(knd) arithmetic
c               cbesf  : vector of ratios of cylindrical Bessel
c                        functions of the first kind with argument
c                        c/2
c               cbesn  : vector of characteristics of these Bessel
c                        functions
c               ibese  : vector of exponents of these Bessel functions
c               cbesdf : vector of ratios of first derivatives of
c                        these Bessel functions
c               cbesdr : vector of ratios of first derivatives of
c                        these Bessel functions to the corresponding
c                        Bessel functions
c               a01    : characteristic of the first a coefficient
c               ia01   : exponent of the first a coefficient
c               b12    : characteristic of the first b coefficient
c               ib12   : exponent of the first b coefficient
c
c     output  : jbpe   : twice the maximum number of terms taken
c               nsub   : maximum number of digits of subtraction
c                        error encountered in the calculation of
c                        the function and its first derivative
c               m1c    : characteristic of Mathieu radial function
c                        of the first kind.
c               im1e   : exponent of Mathieu radial function of the
c                        first kind
c               m1dc   : characteristic of derivative with respect
c                        to x of Mathieu radial function of the first
c                        kind
c               im1de  : exponent of derivative with respect to x of
c                        Mathieu radial function of the first kind
c
        use param
c
c  real(knd) scalars and vectors
        real(knd) a01,b12,cm,coef,dcon,den,dnew,dold,dolds,m1bot,m1c,
     1            m1dc,m1dpos,m1dtemp,m1pos,m1temp,term,termd,termda,
     2            termdb,teste,testeo,u
        real(knd) enr(maxd),enra(maxd),enrb(maxd),cbesdf(maxj),
     1            cbesdr(maxj),cbesf(maxj),cbesn(maxlp)
c
c  integer vectors
        dimension ibese(maxlp)
c
c  convergence ratio dcon is set according to the number of decimal
c  digits available
        dcon=10.0e0_knd**(-ndec-2)
        u=0.5e0_knd*cm
        l2=l/2
c
c  ix=0 for l even, ix=1 for l odd
        ix=l-2*l2
        ics=ix
        if(iopcs.eq.2) ics=iabs(ix-1)
        if(iopcs.eq.1.and.ix.eq.0) iopt=10
        if(iopcs.eq.1.and.ix.eq.1) iopt=11
        if(iopcs.eq.2.and.ix.eq.0) iopt=20
        if(iopcs.eq.2.and.ix.eq.1) iopt=21
        lim=limbpe/2-ix
        nfac=nex-ndec
        teste=10.0e0_knd**nfac
        testeo=1.0e0_knd/teste
        im1dtemp=0
        im1temp=0
        joff=0
        if(iopt.eq.20) joff=1
        jlim=1+joff
        is=l+iss
        if(is.gt.ismax) is=ismax
        lims=is+joff+100+int(cm/10.0)
        limmax=max(lims,lim)
        if(limmax.gt.limd-1) limmax=limd-1
        if(ics.eq.1) go to 20
          do 10 i=1,limmax
          enr(i)=enra(i)
10        continue
          m1bot=a01
          ibot=ia01
          go to 40
20        do 30 i=1,limmax
          enr(i)=enrb(i)
30        continue
          m1bot=b12
          ibot=ib12
40      continue
        den=m1bot
        iden=ibot
          if(is.eq.0) go to 60
          do 50 k=1,is
          den=den*enr(k+joff)
          iterm=int(log10(abs(den)))
          den=den*(10.0e0_knd**(-iterm))
          iden=iden+iterm
50        continue
60      continue
c
c  compute radial function of the first kind m1
c
c    forward summation of product series for m1 and m1d
        coef=cbesn(-l2+is+1+joff)*cbesn(l2+is+1+joff)
        icoef=ibese(-l2+is+1+joff)+ibese(l2+is+1+joff)
        dold=1.0e0_knd
        dolds=dold
        if(iopcs.eq.2) go to 70
        if(iopt.eq.10) term=1.0e0_knd
        if(iopt.eq.11) term=cbesf(l2+is+1)
        m1temp=dold*(term+term)
        m1pos=0.0e0_knd
        if(m1temp.gt.0.0e0_knd) m1pos=m1temp
        go to 80
70      if(iopt.eq.20) termda=cbesdr(is-l2+2)
        if(iopt.eq.20) termdb=-cbesdr(l2+is+2)
        if(iopt.eq.21) termda=cbesdr(is-l2+1)*cbesf(l2+is+1)
        if(iopt.eq.21) termdb=-cbesdr(l2+is+2)*cbesf(l2+is+1)
        termd=dold*(termda+termdb)
        m1dtemp=-u*(termd+termd)
        m1dpos=0.e0_knd
        if(m1dtemp.gt.0.0e0_knd) m1dpos=m1dtemp
80      jbpe=l2
        lflag=0
c
c      forward summation for is => j
          do 100 j=l2+1,is+joff
          dnew=-dold*enr(j)*cbesf(j+is+joff)/cbesf(is+joff-j+1)
          if(iopcs.eq.2) go to 85
          if(iopt.eq.10) term=dnew
          if(iopt.eq.11) term=dnew*cbesf(j+is+1)
          term=term+term
          if(term.gt.0.0e0_knd) m1pos=m1pos+term
          m1temp=m1temp+term
          if(m1temp.eq.0.0e0_knd) go to 120
          if(abs(term/m1temp).lt.dcon) go to 120
          if(abs(m1temp).lt.teste) go to 90
          m1temp=m1temp*testeo
          dnew=dnew*testeo
          dnew=dnew*testeo
          m1pos=m1pos*testeo
          im1temp=im1temp+nfac
          iflag=1
          go to 90
85        if(iopt.eq.20) termda=dnew*cbesdr(is-j+2)
          if(iopt.eq.20) termdb=-dnew*cbesdr(j+is+2)
          if(iopt.eq.21) termda=dnew*cbesdr(is-j+1)*cbesf(j+is+1)
          if(iopt.eq.21) termdb=-dnew*cbesdr(j+is+2)*cbesf(j+is+1)
          termd=termda+termdb
          termd=-u*(termd+termd)
          if(termd.gt.0.0e0_knd) m1dpos=m1dpos+termd
          m1dtemp=m1dtemp+termd
          if(m1dtemp.eq.0.0e0_knd) go to 120
          if(abs(termd/m1dtemp).lt.dcon) go to 120
          if(abs(m1dtemp).lt.teste) go to 90
          m1dtemp=m1dtemp*testeo
          m1dpos=m1dpos*testeo
          im1dtemp=im1dtemp+nfac
          lflag=1
90        dold=dnew
100       continue
c
c      forward summation for j > is
          do 110 j=is+joff+1,lim
          dnew=-dold*enr(j)*cbesf(j-is-joff)*cbesf(j+is+joff)
          if(iopcs.eq.2) go to 105
          if(iopt.eq.10) term=dnew
          if(iopt.eq.11) term=dnew*cbesf(j+is+1)
          term=term+term
          if(term.gt.0.0e0_knd) m1pos=m1pos+term
          m1temp=m1temp+term
          if(m1temp.eq.0.0e0_knd) go to 120
          if(abs(term/m1temp).lt.dcon) go to 120
          go to 110
105       if(iopt.eq.20) termda=dnew*cbesdr(j-is)
          if(iopt.eq.20) termdb=-dnew*cbesdr(j+is+2)
          if(iopt.eq.21) termda=dnew*cbesdr(j-is+1)*cbesf(j+is+1)
          if(iopt.eq.21) termdb=-dnew*cbesdr(j+is+2)*cbesf(j+is+1)
          termd=termda+termdb
          termd=-u*(termd+termd)
          if(termd.gt.0.0e0_knd) m1dpos=m1dpos+termd
          m1dtemp=m1dtemp+termd
          if(m1temp*m1dtemp.eq.0.0e0_knd) go to 120
          if(abs(termd/m1dtemp).lt.dcon) go to 120
110       dold=dnew
120     jbpe=j
c
c  backward summation of product series for m1 and m1d
        if(l2.lt.jlim.or.lflag.eq.1) go to 140
        dold=dolds
          do 130 j=l2,jlim,-1
          dnew=-dold*cbesf(is-j+1+joff)/(enr(j)*cbesf(j+is+joff))
          if(iopcs.eq.2) go to 125
          if(iopt.eq.10) term=dnew
          if(iopt.eq.11) term=dnew*cbesf(j+is)
          term=term+term
          if(term.gt.0.0e0_knd) m1pos=m1pos+term
          m1temp=m1temp+term
          if(m1temp.eq.0.0e0_knd) go to 140
          if(abs(term/m1temp).lt.dcon) go to 140
          go to 130
125       if(iopt.eq.20) termda=dnew*cbesdr(is-j+3)
          if(iopt.eq.20) termdb=-dnew*cbesdr(j+is+1)
          if(iopt.eq.21) termda=dnew*cbesdr(is-j+2)*cbesf(j+is)
          if(iopt.eq.21) termdb=-dnew*cbesdr(j+is+1)*cbesf(j+is)
          termd=termda+termdb
          termd=-u*(termd+termd)
          if(termd.gt.0.0e0_knd) m1dpos=m1dpos+termd
          m1dtemp=m1dtemp+termd
          if(m1dtemp.eq.0.0e0_knd) go to 140
          if(abs(termd/m1dtemp).lt.dcon) go to 140
130       dold=dnew
140     nsub=0
        nsubd=0
        if(iopcs.eq.2) go to 170
        if((m1temp*m1pos).ne.0.0e0_knd) nsub=
     1                             int(log10(abs(m1pos/m1temp)))
        if(nsub.lt.0) nsub=0
        if(nsub.gt.ndec) nsub=ndec
        if(m1temp.eq.0.0e0_knd.and.m1pos.ne.0.0e0_knd) nsub=ndec
        iterm=log10(abs(m1temp))
        m1temp=m1temp*10.0e0_knd**(-iterm)
        m1c=m1temp*coef/den
        if(is.eq.0) m1c=m1c/2.0e0_knd
        im1e=im1temp+iterm+icoef-iden
        iterm=log10(abs(m1c))
        im1e=im1e+iterm
        m1c=m1c*10.0e0_knd**(-iterm)
        if(abs(m1c).ge.1.0e0_knd) go to 150
        m1c=m1c*10.0e0_knd
        im1e=im1e-1
150	continue
        m1c=-m1c
        if(iopt.eq.10) m1c=-m1c
        if(ix.eq.1) m1c=-m1c
        if(2*(is/2).ne.is) m1c=-m1c
c        write(40,160) jbpe,lim,nsub,is
c160     format(3x,'mc1 converged in ',i5,' of ',i5,' terms;',
c     1         ' sub. errors =',i3,'. s = ',i6)
        m1dc=0.0e0_knd
        im1de=0
        go to 200
170     nsubd=0
        if((m1dtemp*m1dpos).ne.0.0e0_knd) nsubd=
     1                             int(log10(abs(m1dpos/m1dtemp)))
        if(nsubd.lt.0) nsubd=0
        if(nsubd.gt.ndec) nsubd=ndec
        if(m1dtemp.eq.0.0e0_knd.and.m1dpos.ne.0.0e0_knd) nsubd=ndec
        iterm=log10(abs(m1dtemp))
        m1dtemp=m1dtemp*10.0e0_knd**(-iterm)
        m1dc=m1dtemp*coef/den
        if(is.eq.0) m1dc=m1dc/2.0e0_knd
        im1de=im1dtemp+iterm+icoef-iden
        iterm=log10(abs(m1dc))
        im1de=im1de+iterm
        m1dc=m1dc*10.0e0_knd**(-iterm)
        if(abs(m1dc).ge.1.0e0_knd) go to 180
        m1dc=m1dc*10.0e0_knd
        im1de=im1de-1
180     continue
        m1dc=-m1dc
        if(iopt.eq.10) m1dc=-m1dc
        if(ix.eq.1) m1dc=-m1dc
        if(2*(is/2).ne.is) m1dc=-m1dc
        m1c=0.0e0_knd
        im1e=0
c        write(40,190) jbpe,lim,nsubd,is
c190     format(3x,'ms1d converged in ',i5,' of ',i5,' terms;',
c     1         ' sub. errors =',i3,'. s = ',i6)
200     if(nsub.eq.0.and.nsubd.eq.0) go to 210
        issinc=(nsub+nsubd)*(1.5+15)+max(0,2*int(log10(cm)))
        if(l.lt.20.and.cm.gt.2000.0e0_knd) issinc=
     1          issinc*cm/2000.0e0_knd
        iss=iss+issinc
210     continue
        return
        end
c
c
        subroutine m2bpe (iopcs,iss,l,c,x1,limbpe,ndec,maxd,enra,
     1                    enrb,maxj,maxlp,maxn,nex,cbesf,cbesn,ibese,
     2                    cbesdf,cbesdr,cneuf,cneun,ineue,cneudf,
     3                    cneudr,m2bot,ibot,jbpe,nsub,m2c,im2e,m2dc,
     4                    im2de)
c
c  purpose    : to calculate the Mathieu radial function of the
c               second kind and its first derivative with respect
c               to the traditional radial coordinate z, using a series
c               of products of cylindrical Bessel functions of the
c               first and second kinds, with the integer Bessel
c               function order offset s (is) for input as a free
c               parameter. c is real here.
c
c  parameters :
c
c     input   : iopcs  : =1 if cosine radial functions calculated
c                      : =2 if sine radial functions calculated
c               iss    : integer offset minus l from calculation
c                        for radial function of the same kind
c                        (sine or cosine) for previous l
c               l      : l
c               c      : c, c being a real positive number
c               x1     : xi-1
c               limbpe : twice maximum number of terms available
c               ndec   : number of decimal digits available in
c                        real(knd) arithmetic
c               maxd   : dimension of enra and enrb vectors
c               enra   : vector of a coefficient ratios
c               enrb   : vector of b coefficient ratios
c               maxj   : dimension of cbesf, cbesdf, cbesdr, cbesf1,
c                        cbesdf1 and cbesdr1 vectors
c               maxlp  : dimension of cbesn, ibese, cbesn1 and ibese1
c                        vectors
c               maxn   : dimension of cneuf, cneudf and cneudr vectors
c               nex    : largest integer exponent abvailable in
c                        real(knd) arithmetic
c               cbesf  : vector of ratios of cylindrical Bessel
c                        functions of the first kind with argument
c                        c*[xi - sqrt(xi*xi - 1)]/2
c               cbesn  : vector of characteristics of these Bessel
c                        functions
c               ibese  : vector of exponents of these Bessel functions
c               cbesdf : vector of ratios of first derivatives of
c                        these Bessel functions
c               cbesdr : vector of ratios of first derivatives of
c                        these Bessel functions to the corresponding
c                        Bessel functions
c               cneuf  : vector of ratios of cylindrical Bessel
c                        functions of the second kind with argument
c                        c*[xi + sqrt(xi*xi - 1)]/2
c               cneun  : vector of characteristics of these Bessel
c                        functions
c               ineue  : vector of exponents of these Bessel functions
c               cneudf : vector of ratios of first derivatives of these
c                        Bessel functions
c               cneudr : vector of ratios of first derivatives of
c                        these Bessel functions to the corresponding
c                        Bessel functions
c               m2bot  : equal to a01, the characteristic of the first
c                        a coefficient when iopcs =1 and equal to b12,
c                        the characteristic of the first b coefficient
c                        when iopcs =2
c               ibot   : equal to ia01, the exponent of the first
c                        a coefficient when iopcs =1 and equal to ib12,
c                        the exponent of the first b coefficient
c                        when iopcs =2
c
c     output  : jbpe   : maximum number of terms taken
c               nsub   : maximum number of digits of subtraction
c                        error encountered in the calculation of
c                        the function and its first derivative
c               m2c    : characteristic of Mathieu radial function
c                        of the second kind.
c               im2e   : exponent of Mathieu radial function of the
c                        second kind
c               m2dc   : characteristic of derivative with respect
c                        to x of Mathieu radial function of the second
c                        kind
c               im2de  : exponent of derivative with respect to x of
c                        Mathieu radial function of the second kind
c
        use param
c
c  real(knd) scalars and vectors
        real(knd) c,coef,coefa,coefb,dcon,den,dnewa,dnewb,dolda,
     1            doldas,doldb,doldbs,m2bot,m2c,m2dc,m2dpos,m2dtemp,
     2            m2pos,m2temp,sqx2m1,term,terma,termb,termd,
     3            termd1,termd1a,termd1b,termd2,termd2a,termd2b,teste,
     4            testeo,u1,u2,xi,x1,x2m1
        real(knd) enr(maxd),enra(maxd),enrb(maxd),cbesdf(maxj),
     1            cbesdr(maxj),cbesf(maxj),cbesn(maxlp),cneudf(maxn),
     2            cneudr(maxn),cneuf(maxn),cneun(maxlp)
c
c  integer vectors
        dimension ibese(maxlp),ineue(maxlp)
c
c  convergence ratio dcon is set according to the number of decimal
c  digits available
        dcon=10.0e0_knd**(-ndec-2)
        xi=x1+1.0e0_knd
        x2m1=x1*(x1+2.0e0_knd)
        sqx2m1=sqrt(x2m1)
        u1=0.5e0_knd*c/(xi+sqx2m1)
        u2=0.5e0_knd*c*(xi+sqx2m1)
        l2=l/2
c
c  ix=0 for l even, ix=1 for l odd
        ix=l-2*l2
        ics=ix
        if(iopcs.eq.2) ics=iabs(ix-1)
        if(iopcs.eq.1.and.ix.eq.0) iopt=10
        if(iopcs.eq.1.and.ix.eq.1) iopt=11
        if(iopcs.eq.2.and.ix.eq.0) iopt=20
        if(iopcs.eq.2.and.ix.eq.1) iopt=21
        lim=limbpe/2-ix
        nexpt=nex-ndec
        nfac=50
        teste=10.0e0_knd**(nex-10)
        testeo=10.0e0_knd**(-nfac)
        im2dtemp=0
        im2temp=0
        if(iopcs.eq.2) go to 30
          do 20 i=1,lim
          enr(i)=enra(i)
20        continue
          go to 50
30        do 40 i=1,lim
          enr(i)=enrb(i)
40        continue
50      continue
        joff=0
        if(iopt.eq.20) joff=1
        jlim=1+joff
        den=m2bot
        iden=ibot
        if(iss.gt.(l2-joff)) iss=l2-joff
        is=iss
        if(is.eq.0) go to 70
          do 60 k=1,is
          den=den*enr(k+joff)
          iterm=int(log10(abs(den)))
          den=den*(10.0e0_knd**(-iterm))
60        iden=iden+iterm
70      continue
c
c  compute radial function of the second kind m2
c
c  forward summation of product series for m2 and m2d
        coefa=cbesn(l2-is+1-joff)*cneun(l2+is+1+joff)
        inda=ibese(l2-is+1-joff)+ineue(l2+is+1+joff)
        coefb=cbesn(l2+is+1+joff)*cneun(l2-is+1-joff)
        indb=ibese(l2+is+1+joff)+ineue(l2-is+1-joff)
        if(inda.lt.indb) go to 80
        coef=coefa
        icoef=inda
        dolda=1.0e0_knd
        indexp=indb-inda
          if(indexp.lt.-nexpt) then
          dolda=dolda*(10.0e0_knd**(nexpt))
          indexp=indexp+nexpt
          icoef=icoef-nexpt
          end if
          if(indexp.lt.-nexpt) then
          doldb=0.0e0_knd
          else
          doldb=(coefb/coefa)*(10.0e0_knd**(indexp))
          end if
        go to 90
80      doldb=1.0e0_knd
        dolda=(coefa/coefb)*(10.0e0_knd**(inda-indb))
        coef=coefb
        icoef=indb
90      if(is.eq.0.and.iopt.eq.10) dolda=0.5e0_knd*dolda
        if(is.eq.0.and.iopt.eq.10) doldb=0.5e0_knd*doldb
        doldas=dolda
        doldbs=doldb
        if(iopt.eq.10) terma=1.0e0_knd
        if(iopt.eq.10) termb=1.0e0_knd
        if(iopt.eq.11) terma=cneuf(l2+is+1)
        if(iopt.eq.11) termb=cbesf(l2+is+1)
        if(iopt.eq.21) terma=cneuf(l2+is+1)
        if(iopt.eq.21) termb=-cbesf(l2+is+1)
        if(iopt.eq.20) terma=1.0e0_knd
        if(iopt.eq.20) termb=-1.0e0_knd
        m2temp=dolda*terma+doldb*termb
        if(iopt.eq.10) termd1a=cbesdr(l2-is+1)
        if(iopt.eq.10) termd2a=cneudr(l2+is+1)
        if(iopt.eq.10) termd1b=cbesdr(l2+is+1)
        if(iopt.eq.10) termd2b=cneudr(l2-is+1)
        if(iopt.eq.20) termd1a=cbesdr(l2-is)
        if(iopt.eq.20) termd2a=cneudr(l2+is+2)
        if(iopt.eq.20) termd1b=-cbesdr(l2+is+2)
        if(iopt.eq.20) termd2b=-cneudr(l2-is)
        if(iopt.eq.11) termd1a=cbesdr(l2-is+1)*cneuf(l2+is+1)
        if(iopt.eq.11) termd2a=cneudr(l2+is+2)*cneuf(l2+is+1)
        if(iopt.eq.11) termd1b=cbesdr(l2+is+2)*cbesf(l2+is+1)
        if(iopt.eq.11) termd2b=cneudr(l2-is+1)*cbesf(l2+is+1)
        if(iopt.eq.21) termd1a=cbesdr(l2-is+1)*cneuf(l2+is+1)
        if(iopt.eq.21) termd2a=cneudr(l2+is+2)*cneuf(l2+is+1)
        if(iopt.eq.21) termd1b=-cbesdr(l2+is+2)*cbesf(l2+is+1)
        if(iopt.eq.21) termd2b=-cneudr(l2-is+1)*cbesf(l2+is+1)
        termd1=dolda*termd1a+doldb*termd1b
        termd2=dolda*termd2a+doldb*termd2b
        m2dtemp=-u1*termd1+u2*termd2
        m2pos=0.0e0_knd
        if(m2temp.gt.0.0e0_knd) m2pos=m2temp
        m2dpos=0.e0_knd
        if(m2dtemp.gt.0.0e0_knd) m2dpos=m2dtemp
        jbpe=l2
          do 100 j=l2+1,lim
          dnewa=-dolda*enr(j)*cbesf(j-is-joff)*cneuf(j+is+joff)
          dnewb=-doldb*enr(j)*cbesf(j+is+joff)*cneuf(j-is-joff)
          if(iopt.eq.10) terma=dnewa
          if(iopt.eq.10) termb=dnewb
          if(iopt.eq.11) terma=dnewa*cneuf(j+is+1)
          if(iopt.eq.11) termb=dnewb*cbesf(j+is+1)
          if(iopt.eq.21) terma=dnewa*cneuf(j+is+1)
          if(iopt.eq.21) termb=-dnewb*cbesf(j+is+1)
          if(iopt.eq.20) terma=dnewa
          if(iopt.eq.20) termb=-dnewb
          term=terma+termb
          if(iopt.eq.10) termd1a=dnewa*cbesdr(j-is+1)
          if(iopt.eq.10) termd2a=dnewa*cneudr(j+is+1)
          if(iopt.eq.10) termd1b=dnewb*cbesdr(j+is+1)
          if(iopt.eq.10) termd2b=dnewb*cneudr(j-is+1)
          if(iopt.eq.20) termd1a=dnewa*cbesdr(j-is)
          if(iopt.eq.20) termd2a=dnewa*cneudr(j+is+2)
          if(iopt.eq.20) termd1b=-dnewb*cbesdr(j+is+2)
          if(iopt.eq.20) termd2b=-dnewb*cneudr(j-is)
          if(iopt.eq.11) termd1a=dnewa*cbesdr(j-is+1)*cneuf(j+is+1)
          if(iopt.eq.11) termd2a=dnewa*cneudr(j+is+2)*cneuf(j+is+1)
          if(iopt.eq.11) termd1b=dnewb*cbesdr(j+is+2)*cbesf(j+is+1)
          if(iopt.eq.11) termd2b=dnewb*cneudr(j-is+1)*cbesf(j+is+1)
          if(iopt.eq.21) termd1a=dnewa*cbesdr(j-is+1)*cneuf(j+is+1)
          if(iopt.eq.21) termd2a=dnewa*cneudr(j+is+2)*cneuf(j+is+1)
          if(iopt.eq.21) termd1b=-dnewb*cbesdr(j+is+2)*cbesf(j+is+1)
          if(iopt.eq.21) termd2b=-dnewb*cneudr(j-is+1)*cbesf(j+is+1)
          termd1=termd1a+termd1b
          termd2=termd2a+termd2b
          termd=-u1*termd1+u2*termd2
          if(term.gt.0.0e0_knd) m2pos=m2pos+term
          if(termd.gt.0.0e0_knd) m2dpos=m2dpos+termd
          m2temp=m2temp+term
          m2dtemp=m2dtemp+termd
          if(m2temp*m2dtemp.eq.0.0e0_knd) go to 110
          if((abs(term/m2temp)+abs(termd/m2dtemp)).lt.dcon) go to 110
            if(abs(m2temp).gt.teste.or.abs(m2dtemp).gt.teste) then
            m2temp=m2temp*testeo
            m2dtemp=m2dtemp*testeo
            dnewa=dnewa*testeo
            dnewb=dnewb*testeo
            m2pos=m2pos*testeo
            m2dpos=m2dpos*testeo
            doldas=doldas*testeo
            doldbs=doldbs*testeo
            im2temp=im2temp+nfac
            im2dtemp=im2dtemp+nfac
            end if
          dolda=dnewa
          doldb=dnewb
100       continue
110     jbpe=j
c
c  backward summation of m2 product series (terms with j > is)
        if (l2.eq.0) go to 170
        dolda=doldas
        doldb=doldbs
        jlim1=jlim+is
        if(jlim1.gt.l2) go to 140
        do 130 j=l2,jlim1,-1
          dnewa=-dolda/(enr(j)*cbesf(j-is-joff)*cneuf(j+is+joff))
          dnewb=-doldb/(enr(j)*cbesf(j+is+joff)*cneuf(j-is-joff))
          if(iopt.eq.10) terma=dnewa
          if(iopt.eq.10) termb=dnewb
          if(iopt.eq.11) terma=dnewa*cneuf(j+is)
          if(iopt.eq.11) termb=dnewb*cbesf(j+is)
          if(iopt.eq.21) terma=dnewa*cneuf(j+is)
          if(iopt.eq.21) termb=-dnewb*cbesf(j+is)
          if(iopt.eq.20) terma=dnewa
          if(iopt.eq.20) termb=-dnewb
          term=terma+termb
          if(iopt.eq.10) termd1a=dnewa*cbesdr(j-is)
          if(iopt.eq.10) termd2a=dnewa*cneudr(j+is)
          if(iopt.eq.10) termd1b=dnewb*cbesdr(j+is)
          if(iopt.eq.10) termd2b=dnewb*cneudr(j-is)
          if(iopt.eq.20) termd1a=dnewa*cbesdr(j-is-1)
          if(iopt.eq.20) termd2a=dnewa*cneudr(j+is+1)
          if(iopt.eq.20) termd1b=-dnewb*cbesdr(j+is+1)
          if(iopt.eq.20) termd2b=-dnewb*cneudr(j-is-1)
          if(iopt.eq.11) termd1a=dnewa*cbesdr(j-is)*cneuf(j+is)
          if(iopt.eq.11) termd2a=dnewa*cneudr(j+is+1)*cneuf(j+is)
          if(iopt.eq.11) termd1b=dnewb*cbesdr(j+is+1)*cbesf(j+is)
          if(iopt.eq.11) termd2b=dnewb*cneudr(j-is)*cbesf(j+is)
          if(iopt.eq.21) termd1a=dnewa*cbesdr(j-is)*cneuf(j+is)
          if(iopt.eq.21) termd2a=dnewa*cneudr(j+is+1)*cneuf(j+is)
          if(iopt.eq.21) termd1b=-dnewb*cbesdr(j+is+1)*cbesf(j+is)
          if(iopt.eq.21) termd2b=-dnewb*cneudr(j-is)*cbesf(j+is)
          termd1=termd1a+termd1b
          termd2=termd2a+termd2b
          termd=-u1*termd1+u2*termd2
          if(term.gt.0.0e0_knd) m2pos=m2pos+term
          if(termd.gt.0.0e0_knd) m2dpos=m2dpos+termd
          m2temp=m2temp+term
          m2dtemp=m2dtemp+termd
          if(m2temp*m2dtemp.eq.0.0e0_knd) go to 140
          if((abs(term/m2temp)+abs(termd/m2dtemp)).lt.dcon) go to 140
            if(abs(m2temp).gt.teste.or.abs(m2dtemp).gt.teste) then
            m2temp=m2temp*testeo
            m2dtemp=m2dtemp*testeo
            dnewa=dnewa*testeo
            dnewb=dnewb*testeo
            m2pos=m2pos*testeo
            m2dpos=m2dpos*testeo
            im2temp=im2temp+nfac
            im2dtemp=im2dtemp+nfac
            end if
120       dolda=dnewa
          doldb=dnewb
130       continue
140     continue
c  backward summation of m2 product series (terms with j <= is)
        if(is.eq.0) go to 170
        jlow=jlim1-1
        if(jlim1.gt.l2) jlow=l2
          do 160 j=jlow,jlim,-1
          dnewa=dolda*cbesf(is-j+1+joff)/(enr(j)*cneuf(j+is+joff))
          dnewb=doldb*cneuf(is-j+1+joff)/(enr(j)*cbesf(j+is+joff))
          if(iopt.eq.10) terma=dnewa
          if(iopt.eq.10) termb=dnewb
          if(iopt.eq.11) terma=dnewa*cneuf(j+is)
          if(iopt.eq.11) termb=dnewb*cbesf(j+is)
          if(iopt.eq.21) terma=dnewa*cneuf(j+is)
          if(iopt.eq.21) termb=-dnewb*cbesf(j+is)
          if(iopt.eq.20) terma=dnewa
          if(iopt.eq.20) termb=-dnewb
          term=terma+termb
          if(iopt.eq.10) termd1a=dnewa*cbesdr(is-j+2)
          if(iopt.eq.10) termd2a=dnewa*cneudr(j+is)
          if(iopt.eq.10) termd1b=dnewb*cbesdr(j+is)
          if(iopt.eq.10) termd2b=dnewb*cneudr(is-j+2)
          if(iopt.eq.20) termd1a=dnewa*cbesdr(is-j+3)
          if(iopt.eq.20) termd2a=dnewa*cneudr(j+is+1)
          if(iopt.eq.20) termd1b=-dnewb*cbesdr(j+is+1)
          if(iopt.eq.20) termd2b=-dnewb*cneudr(is-j+3)
          if(iopt.eq.11) termd1a=dnewa*cbesdr(is-j+2)*cneuf(j+is)
          if(iopt.eq.11) termd2a=dnewa*cneudr(j+is+1)*cneuf(j+is)
          if(iopt.eq.11) termd1b=dnewb*cbesdr(j+is+1)*cbesf(j+is)
          if(iopt.eq.11) termd2b=dnewb*cneudr(is-j+2)*cbesf(j+is)
          if(iopt.eq.21) termd1a=dnewa*cbesdr(is-j+2)*cneuf(j+is)
          if(iopt.eq.21) termd2a=dnewa*cneudr(j+is+1)*cneuf(j+is)
          if(iopt.eq.21) termd1b=-dnewb*cbesdr(j+is+1)*cbesf(j+is)
          if(iopt.eq.21) termd2b=-dnewb*cneudr(is-j+2)*cbesf(j+is)
          termd1=termd1a+termd1b
          termd2=termd2a+termd2b
          termd=-u1*termd1+u2*termd2
          if(term.gt.0.0e0_knd) m2pos=m2pos+term
          if(termd.gt.0.0e0_knd) m2dpos=m2dpos+termd
          m2temp=m2temp+term
          m2dtemp=m2dtemp+termd
          if(m2temp*m2dtemp.eq.0.0e0_knd) go to 170
          if((abs(term/m2temp)+abs(termd/m2dtemp)).lt.dcon) go to 170
            if(abs(m2temp).gt.teste.or.abs(m2dtemp).gt.teste) then
            m2temp=m2temp*testeo
            m2dtemp=m2dtemp*testeo
            dnewa=dnewa*testeo
            dnewb=dnewb*testeo
            m2pos=m2pos*testeo
            m2dpos=m2dpos*testeo
            im2temp=im2temp+nfac
            im2dtemp=im2dtemp+nfac
            end if
150       dolda=dnewa
          doldb=dnewb
160       continue
170     nsub=0
        if((m2temp*m2pos).ne.0.0e0_knd) nsub=
     1                             int(log10(abs(m2pos/m2temp)))
        if(nsub.lt.0) nsub=0
        if(nsub.gt.ndec) nsub=ndec
        if(m2temp.eq.0.0e0_knd.and.m2pos.ne.0.0e0_knd) nsub=ndec
        nsubd=0
        if((m2dtemp*m2dpos).ne.0.0e0_knd) nsubd=
     1                             int(log10(abs(m2dpos/m2dtemp)))
        if(nsubd.lt.0) nsubd=0
        if(nsubd.gt.ndec) nsubd=ndec
        if(m2dtemp.eq.0.0e0_knd.and.m2dpos.ne.0.0e0_knd) nsubd=ndec
c        if(iopcs.eq.1) write(40,180) jbpe,lim,nsub,nsubd,is
c180     format(3x,'mc2/mc2d converged in ',i5,' of ',i5,' terms;',
c     1         ' sub. errors =',i3,' and',i3,' digits. s = ',i6)
c        if(iopcs.eq.2) write(40,190) jbpe,lim,nsub,nsubd,is
c190     format(3x,'ms2/ms2d converged in ',i5,' of ',i5,' terms;',
c     1         ' sub. errors =',i3,' and',i3,' digits. s = ',i6)
        iss=is
        llim1=int(0.6e0_knd*c)
        llim2=int(2.0e0_knd*c/3.1416)
        if((nsub.gt.2.and.nsubd.gt.2).and.l.ge.llim1)
     1                 iss=is+max(2,int(abs(log10(c))))
        if((nsub.gt.0.and.nsubd.gt.0).and.l.ge.llim2)
     1                 iss=is+max(2,int(abs(log10(c))))
        iterm=log10(abs(m2temp))
        m2temp=m2temp*10.0e0_knd**(-iterm)
        m2c=m2temp*coef/den
        im2e=im2temp+iterm+icoef-iden
        iterm=log10(abs(m2c))
        im2e=im2e+iterm
        m2c=m2c*10.0e0_knd**(-iterm)
        if(abs(m2c).ge.1.0e0_knd) go to 200
        m2c=m2c*10.0e0_knd
        im2e=im2e-1
200	continue
        iterm=log10(abs(m2dtemp))
        m2dtemp=m2dtemp*10.0e0_knd**(-iterm)
        m2dc=m2dtemp*coef/den
        im2de=im2dtemp+iterm+icoef-iden
        iterm=log10(abs(m2dc))
        im2de=im2de+iterm
        m2dc=m2dc*10.0e0_knd**(-iterm)
        if(abs(m2dc).ge.1.0e0_knd) go to 210
        m2dc=m2dc*10.0e0_knd
        im2de=im2de-1
210     continue
        return
        end
c
c
        subroutine m3bpe (iopcs,l,cm,x1,limbpe,ndec,maxd,enra,enrb,
     1                    maxj,maxlp,maxn,nex,cbesf,cbesn,ibese,
     2                    cbesdf,cbesdr,cneuf,cneun,ineue,cneudf,
     3                    cneudr,a01,ia01,b12,ib12,pi,jbpe,nsub,
     4                    m3c,im3e,m3dc,im3de)
c
c  purpose    : to calculate the Mathieu radial function of the
c               third kind and its first derivative with respect
c               to the traditional radial coordinate z, for the
c               case where q is negative, using a series of products
c               of modified cylindrical Bessel functions of the first
c               and second kinds
c
c  parameters :
c
c     input   : iopcs  : =1 if cosine radial functions calculated
c                      : =2 if sine radial functions calculated
c               l      : l
c               cm     : magnitude of c (c being positive imaginary)
c               x1     : xi-1
c               limbpe : twice maximum number of terms available
c               ndec   : number of decimal digits available in
c                        real(knd) arithmetic
c               maxd   : dimension of enra and enrb vectors
c               enra   : vector of a coefficient ratios
c               enrb   : vector of b coefficient ratios
c               maxj   : dimension of cbesf, cbesdf and cbesdr vectors
c               maxlp  : dimension of cbesn, ibese, cneun and ineue
c                        vectors
c               maxn   : dimension of cneuf, cneudf and cneudr vectors
c               nex    : largest integer exponent abvailable in
c                        real(knd) arithmetic
c               cbesf  : vector of ratios of modified cylindrical Bessel
c                        functions of the first kind with argument
c                        c*[xi - sqrt(xi*xi - 1)]/2
c               cbesn  : vector of characteristics of these Bessel
c                        functions
c               ibese  : vector of exponents of these Bessel functions
c               cbesdf : vector of ratios of first derivatives of
c                        these Bessel functions
c               cbesdr : vector of ratios of first derivatives of
c                        these Bessel functions to the corresponding
c                        Bessel functions
c               cneuf  : vector of ratios of modified cylindrical Bessel
c                        functions of the second kind with argument
c                        c*[xi + sqrt(xi*xi - 1)]/2
c               cneun  : vector of characteristics of these Bessel
c                        functions
c               ineue  : vector of exponents of these Bessel functions
c               cneudf : vector of ratios of first derivatives of these
c                        Bessel functions
c               cneudr : vector of ratios of first derivatives of
c                        these Bessel functions to the corresponding
c                        Bessel functions
c               a01    : characteristic of the first a coefficient
c               ia01   : exponent of the first a coefficient
c               b12    : characteristic of the first b coefficient
c               ib12   : exponent of the first b coefficient
c               pi     : 3.14...
c
c     output  : jbpe   : maximum number of terms taken
c               nsub   : maximum number of digits of subtraction
c                        error encountered in the calculation of
c                        the function and its first derivative
c               m3c    : characteristic of Mathieu radial function
c                        of the third kind [note that when l is even,
c                        the Mathieu radial function of the third kind
c                        is equal to i times m3c times 10**(im3e)]
c               im3e   : exponent of Mathieu radial function of the
c                        third kind
c               m3dc   : characteristic of derivative with respect
c                        to x of Mathieu radial function of the third
c                        kind (see note above for m3c)
c               im3de  : exponent of derivative with respect to x of
c                        Mathieu radial function of the third kind
c
        use param
c
c  real(knd) scalars and vectors
        real(knd) a01,b12,cm,dcon,dnew,dold,m3bot,m3c,m3dc,
     1            m3dtemp,m3dpos,m3pos,m3temp,pi,sqx2m1,term,termd,
     2            termd1,termd2,teste,testeo,u1,u2,xi,x1,x2m1
        real(knd) enr(maxd),enra(maxd),enrb(maxd),cbesdf(maxj),
     1            cbesdr(maxj),cbesf(maxj),cbesn(maxlp),cneudf(maxn),
     2            cneudr(maxn),cneuf(maxn),cneun(maxlp)
c
c  integer vectors
        dimension ibese(maxlp),ineue(maxlp)
c
c  convergence ratio dcon is set according to the number of decimal
c  digits available
        dcon=10.0e0_knd**(-ndec-2)
        xi=x1+1.0e0_knd
        x2m1=x1*(x1+2.0e0_knd)
        sqx2m1=sqrt(x2m1)
        u1=0.5e0_knd*cm/(xi+sqx2m1)
        u2=0.5e0_knd*cm*(xi+sqx2m1)
        l2=l/2
c
c  ix=0 for l even, ix=1 for l odd
        ix=l-2*l2
c  lx2=0 for l2 even, lx2=1 for l2 odd
        lx2=l2-2*(l2/2)
        ics=ix
        if(iopcs.eq.2) ics=iabs(ix-1)
        if(iopcs.eq.1.and.ix.eq.0) iopt=10
        if(iopcs.eq.1.and.ix.eq.1) iopt=11
        if(iopcs.eq.2.and.ix.eq.0) iopt=20
        if(iopcs.eq.2.and.ix.eq.1) iopt=21
        lim=limbpe/2-ix
        nfac=nex-ndec
        teste=10.0e0_knd**nfac
        testeo=1.0e0_knd/teste
        im3dtemp=0
        im3temp=0
        if(ics.eq.1) go to 20
          do 10 i=1,lim
          enr(i)=enra(i)
10        continue
          m3bot=a01
          ibot=ia01
          go to 40
20        do 30 i=1,lim
          enr(i)=enrb(i)
30        continue
          m3bot=b12
          ibot=ib12
40      continue
c
c
c  compute radial function of the third kind m3
c
c  forward summation of product series for m3 and m3d
        m3pos=0.0e0_knd
        m3dpos=0.0e0_knd
        dold=1.0e0_knd
        m3temp=1.0e0_knd
        if(iopt.eq.11) m3temp=cneuf(l2+1)-cbesf(l2+1)
        if(iopt.eq.21) m3temp=cneuf(l2+1)+cbesf(l2+1)
        if(iopt.eq.20) m3temp=cbesf(l2+1)/cneuf(l2)-
     1                        cneuf(l2+1)/cbesf(l2)
        if(iopt.eq.10) termd1=cbesdr(l2+1)
        if(iopt.eq.11) termd1=cbesdr(l2+1)*cneuf(l2+1)-
     1                        cbesdr(l2+2)*cbesf(l2+1)
        if(iopt.eq.21) termd1=cbesdr(l2+1)*cneuf(l2+1)+
     1                        cbesdr(l2+2)*cbesf(l2+1)
        if(iopt.eq.20) termd1=-cbesdr(l2)*cneuf(l2+1)/
     1                        cbesf(l2)+cbesdr(l2+2)*
     2                        cbesf(l2+1)/cneuf(l2)
        if(iopt.eq.10) termd2=cneudr(l2+1)
        if(iopt.eq.11) termd2=cneudr(l2+2)*cneuf(l2+1)-
     1                        cneudr(l2+1)*cbesf(l2+1)
        if(iopt.eq.21) termd2=cneudr(l2+2)*cneuf(l2+1)+
     1                        cneudr(l2+1)*cbesf(l2+1)
        if(iopt.eq.20) termd2=-cneudr(l2+2)*cneuf(l2+1)/
     1                        cbesf(l2)+cneudr(l2)*
     2                        cbesf(l2+1)/cneuf(l2)
        m3dtemp=-u1*termd1+u2*termd2
        if(m3temp.gt.0.0e0_knd) m3pos=m3temp
        if(m3dtemp.gt.0.0e0_knd) m3dpos=m3dtemp
        jbpe=l2
          do 50 j=l2+1,lim
          dnew=dold*enr(j)*cbesf(j)*cneuf(j)
          if(iopt.eq.10) term=dnew
          if(iopt.eq.11) term=dnew*(-cbesf(j+1)+cneuf(j+1))
          if(iopt.eq.21) term=dnew*(cbesf(j+1)+cneuf(j+1))
          if(iopt.eq.20) term=dnew*(-cneuf(j+1)/cbesf(j)+
     1                        cbesf(j+1)/cneuf(j))
          if(iopt.eq.10) termd1=dnew*cbesdr(j+1)
          if(iopt.eq.11) termd1=dnew*(cbesdr(j+1)*cneuf(j+1)
     1                          -cbesdr(j+2)*cbesf(j+1))
          if(iopt.eq.21) termd1=dnew*(cbesdr(j+1)*cneuf(j+1)
     1                          +cbesdr(j+2)*cbesf(j+1))
          if(iopt.eq.20) termd1=dnew*(-cbesdr(j)*cneuf(j+1)/
     1                          cbesf(j)+cbesdr(j+2)*
     2                          cbesf(j+1)/cneuf(j))
          if(iopt.eq.10) termd2=dnew*cneudr(j+1)
          if(iopt.eq.11) termd2=dnew*(cneudr(j+2)*cneuf(j+1)
     1                          -cneudr(j+1)*cbesf(j+1))
          if(iopt.eq.21) termd2=dnew*(cneudr(j+2)*cneuf(j+1)
     1                          +cneudr(j+1)*cbesf(j+1))
          if(iopt.eq.20) termd2=dnew*(-cneudr(j+2)*
     1                          cneuf(j+1)/cbesf(j)+
     2                          cneudr(j)*cbesf(j+1)/
     3                          cneuf(j))
          termd=-u1*termd1+u2*termd2
          if(term.gt.0.0e0_knd) m3pos=m3pos+term
          if(termd.gt.0.0e0_knd) m3dpos=m3dpos+termd
          m3temp=m3temp+term
          m3dtemp=m3dtemp+termd
          if(m3temp*m3dtemp.eq.0.0e0_knd) go to 60
          if((abs(term/m3temp)+abs(termd/m3dtemp)).lt.dcon) go to 60
          dold=dnew
50        continue
60      jbpe=j
c
c  backward summation of product series for m3
        if (l2.eq.0) go to 90
        dold=1.0e0_knd
        jlim=1
        if(iopcs.eq.2.and.ix.eq.0) jlim=2
        if(l2.lt.jlim) go to 90
          do 80 j=l2,jlim,-1
          dnew=dold/(enr(j)*cbesf(j)*cneuf(j))
          if(iopt.eq.10) term=dnew
          if(iopt.eq.11) term=dnew*(-cbesf(j)+cneuf(j))
          if(iopt.eq.21) term=dnew*(cbesf(j)+cneuf(j))
          if(iopt.eq.20) term=dnew*(-cneuf(j)/cbesf(j-1)+
     1                        cbesf(j)/cneuf(j-1))
          if(iopt.eq.10) termd1=dnew*cbesdr(j)
          if(iopt.eq.11) termd1=dnew*(cbesdr(j)*cneuf(j)
     1                          -cbesdr(j+1)*cbesf(j))
          if(iopt.eq.21) termd1=dnew*(cbesdr(j)*cneuf(j)
     1                          +cbesdr(j+1)*cbesf(j))
          if(iopt.eq.20) termd1=dnew*(-cbesdr(j-1)*cneuf(j)/
     1                          cbesf(j-1)+cbesdr(j+1)*
     2                          cbesf(j)/cneuf(j-1))
          if(iopt.eq.10) termd2=dnew*cneudr(j)
          if(iopt.eq.11) termd2=dnew*(cneudr(j+1)*cneuf(j)
     1                          -cneudr(j)*cbesf(j))
          if(iopt.eq.21) termd2=dnew*(cneudr(j+1)*cneuf(j)
     1                          +cneudr(j)*cbesf(j))
          if(iopt.eq.20) termd2=dnew*(-cneudr(j+1)*cneuf(j)/
     1                          cbesf(j-1)+cneudr(j-1)*
     2                          cbesf(j)/cneuf(j-1))
          termd=-u1*termd1+u2*termd2
          if(term.gt.0.0e0_knd) m3pos=m3pos+term
          if(termd.gt.0.0e0_knd) m3dpos=m3dpos+termd
          m3temp=m3temp+term
          m3dtemp=m3dtemp+termd
          if(m3temp*m3dtemp.eq.0.0e0_knd) go to 90
          if((abs(term/m3temp)+abs(termd/m3dtemp)).lt.dcon) go to 90
            if(abs(m3temp).gt.teste.or.abs(m3dtemp).gt.teste) then
            m3temp=m3temp*testeo
            m3dtemp=m3dtemp*testeo
            dnew=dnew*testeo
            m3pos=m3pos*testeo
            m3dpos=m3dpos*testeo
            im3temp=im3temp+nfac
            im3dtemp=im3dtemp+nfac
            end if
70        dold=dnew
80        continue
90      nsub=0
        if((m3temp*m3pos).ne.0.0e0_knd) nsub=
     1                             int(log10(abs(m3pos/m3temp)))
        if(nsub.lt.0) nsub=0
        if(nsub.gt.ndec) nsub=ndec
        if(m3temp.eq.0.0e0_knd.and.m3pos.ne.0.0e0_knd) nsub=ndec
        nsubd=0
        if((m3dtemp*m3dpos).ne.0.0e0_knd) nsubd=
     1                             int(log10(abs(m3dpos/m3dtemp)))
        if(nsubd.lt.0) nsubd=0
        if(nsubd.gt.ndec) nsubd=ndec
        if(m3dtemp.eq.0.0e0_knd.and.m3dpos.ne.0.0e0_knd) nsubd=ndec
        nsub=max(nsub,nsubd)
c        if(iopcs.eq.1) write(40,100) jbpe,lim,nsub
c100     format(3x,'mc3 numerator converged in ',i5,' terms; ',i5,
c     1         ' available; sub. error = ',i3,' digits.')
c        if(iopcs.eq.2) write(40,110) jbpe,lim,nsub
c110     format(3x,'ms3 numerator converged in ',i5,' terms; ',i5,
c     1         ' available; sub. error = ',i3,' digits.')
        iterm=log10(abs(m3temp))
        m3temp=m3temp*10.0e0_knd**(-iterm)
        m3c=2.0e0_knd*m3temp*cbesn(l2+1)*cneun(l2+1)/(m3bot*pi)
        im3e=im3temp+iterm+ibese(l2+1)+ineue(l2+1)-ibot
        iterm=0
        if(m3c.ne.0.0e0_knd) iterm=int(log10(abs(m3c)))
        im3e=im3e+iterm
        m3c=m3c*10.0e0_knd**(-iterm)
        if(abs(m3c).ge.1.0e0_knd) go to 120
        m3c=m3c*10.0e0_knd
        im3e=im3e-1
120	continue
        iterm=log10(abs(m3dtemp))
        m3dtemp=m3dtemp*10.0e0_knd**(-iterm)
        m3dc=2.0e0_knd*m3dtemp*cbesn(l2+1)*cneun(l2+1)/(m3bot*pi)
        im3de=im3dtemp+iterm+ibese(l2+1)+ineue(l2+1)-ibot
        iterm=0
        if(m3dc.ne.0.0e0_knd) iterm=int(log10(abs(m3dc)))
        im3de=im3de+iterm
        m3dc=m3dc*10.0e0_knd**(-iterm)
        if(abs(m3dc).ge.1.0e0_knd) go to 130
        m3dc=m3dc*10.0e0_knd
        im3de=im3de-1
130     continue
        if(lx2.eq.0) m3c=-m3c
        if(lx2.eq.0) m3dc=-m3dc
        if(iopt.eq.20) m3c=-m3c
        if(iopt.eq.20) m3dc=-m3dc
        return
        end
c
c
        subroutine m3neu (iopcs,l,cm,x1,limneu,ndec,minacc,maxd,
     1                    enra,enrb,maxn,maxlp,nex,cneuf,cneun,ineue,
     2                    cneudf,cneudr,m3bot,pi,jneu,ndig,m3c,
     3                    im3e,m3dc,im3de,jtermflag)
c
c  purpose    : to calculate the Mathieu radial function of the
c               third kind (q negative) and its first derivative
c               with respect to the traditional radial coordinate
c               z, using an expansion of modified cylindrical Bessel
c               functions of the second kind with argument c*xi.
c
c  parameters :
c
c     input   : iopcs  : =1 if cosine radial functions calculated
c                      : =2 if sine radial functions calculated
c               l      : l
c               cm     : magnitude of c (c being positive imaginary)
c               x1     : xi-1
c               limneu : twice maximum number of terms available
c               ndec   : number of decimal digits available in
c                        real(knd) arithmetic
c               minacc : desired minimum accuracy of the calculated
c                        functions, in decimal digits
c               maxd   : size of enr vector
c               enra   : a coefficient ratios
c               enrb   : b coefficient ratios
c               maxn   : size of cneuf, cneudf and cneudr vectors
c               maxlp  : maximum l value desired; size of the cneun and
c                        ineue vectors
c               nex    : largest integer exponent abvailable in
c                        real(knd) arithmetic
c               cneuf  : ratios of cylindrical Bessel functions of the
c                        second kind
c               cneun  : vector of characteristics of the Bessel
c                        functions
c               ineue  : vector of exponents of the Bessel functions
c               cneudf : ratios of first derivatives of the cylindrical
c                        Bessel functions
c               cneudr : ratios of first derivatives of the cylindrical
c                        Bessel functions to the corresponding Bessel
c                        functions
c               m3bot  : denominator series sum. when q is positive, it
c                        is calculated in dnorma as ce0 when
c                        iopcs = 1 or in dnormb as sed0 when
c                        iopcs = 2. when q is negative it is calculated
c                        as the denominator series in the calculation
c                        of m1 in subroutine mcs1bes
c               pi     : pi
c
c     output  : jneu   : maximum number of terms taken
c               ndig   : minimum number of digits of convergence of
c                        series
c               m3c    : characteristic of Mathieu radial function
c                        of the second kind.
c               im3e   : exponent of Mathieu radial function of the
c                        second kind [note that when q is negative
c                        (isq=-1) the modified Mathieu function m2 =
c                        i times m1c times 10**(im1e), where i is the
c                        square root of -1.]
c               m3dc   : characteristic of derivative with respect
c                        to x of Mathieu radial function of the second
c                        kind
c               im3de  : exponent of derivative with respect to x of
c                        Mathieu radial function of the second kind
c                        [note that when q is negative (isq=-1) the
c                        modified Mathieu function derivative m2d =
c                        i times m1dc times 10**(im1de), where i is the
c                        square root of -1.]
c
        use param
c
c  real(knd) scalars and vectors
        real(knd) cm,con,dcon,dk,dnew,dnewd,dold,doldd,m3bot,m3c,m3dc,
     1            m3dtemp,m3temp,m3top,pi,sqx2m1,sscale,stest,term,
     2            termd,xi,x1,x2m1
        real(knd) enr(maxd),enra(maxd),enrb(maxd),cneudf(maxn),
     1            cneudr(maxn),cneuf(maxn),cneun(maxlp)
c
c  integer vector
        dimension ineue(maxlp)
c
c  convergence ratio dcon is set according to the desired accuracy
c  and the approximate number of terms expected to be taken in the
c  series
        ncon=max(1,-int(log10(x1)))
        dcon=10.0e0_knd**(-minacc-1-ncon)
        xi=x1+1.0e0_knd
        x2m1=x1*(x1+2.0e0_knd)
        sqx2m1=sqrt(x2m1)
        nex2=nex/2
        stest=10.0e0_knd**nex2
        sscale=10.0e0_knd**(-nex2)
        iscale=0
        l2=l/2
c
c  ix=0 for l even, ix=1 for l odd
        ix=l-2*l2
        ics=ix
        if(iopcs.eq.2) ics=iabs(ix-1)
        jlim=1
        if(ix.eq.0.and.iopcs.eq.2) jlim=2
        lim=limneu/2-ix-1
        con=sqx2m1/xi
        if(ics.eq.1) go to 20
          do 10 i=1,lim
          enr(i)=enra(i)
10        continue
          go to 40
20        do 30 i=1,lim
30        enr(i)=enrb(i)
40      continue
c
c  compute radial functions of the third kind m3 and their
c  first derivatives m3d
c
c  forward summation of numerator series for m3 and m3d
        dold=1.0e0_knd
        doldd=cm
        m3temp=dold
        m3dtemp=doldd
        if(iopcs.eq.1) go to 50
        m3temp=real(l,knd)*dold
        m3dtemp=real(l,knd)*doldd
50      jtermflag=0
        kount=1
        kountm=1
          do 90 j=l2+1,lim
          jj=j+j+ix
          dnew=-dold*enr(j)*cneuf(jj)
            if(dnew*dold.gt.0.0e0_knd) then
            kount=kount+1
            else
            kountm=max(kount,kountm)
            kount=1
            end if
          dnewd=-doldd*enr(j)*cneudf(jj)
          term=dnew
          termd=dnewd
          if(iopcs.eq.1) go to 60
          term=dnew*real(jj,knd)
          termd=dnewd*real(jj,knd)
60        m3temp=m3temp+term
          m3dtemp=m3dtemp+termd
          if((abs(term/m3temp)+abs(termd/m3dtemp)).lt.dcon)
     1          go to 100
          if(j.eq.lim) go to 100
          dold=dnew
          doldd=dnewd
          if(abs(m3dtemp).lt.stest) go to 90
          iscale=iscale+nex2
          dold=dold*sscale
          doldd=doldd*sscale
          m3temp=m3temp*sscale
          m3dtemp=m3dtemp*sscale
90        continue
        jtermflag=1
100     jneu=min(j,lim)
        dk=real(kountm,knd)
        ndig=-int(log10(dk*abs(term/m3temp)))
        ndigd=-int(log10(dk*abs(termd/m3dtemp)))
        ndig=min(ndig,ndigd,ndec)
        if(ndig.lt.0) ndig=0
c
c  backward summation of numerator series for m3
        if(iscale.ne.0) go to 130
        if (l2.lt.jlim) go to 130
        dold=1.0e0_knd
        doldd=cm
          do 120 j=l2,jlim,-1
          jj=j+j+ix
          dnew=-dold/(enr(j)*cneuf(jj))
          dnewd=-doldd/(enr(j)*cneudf(jj))
          term=dnew
          termd=dnewd
          if(iopcs.eq.1) go to 110
          term=dnew*real(jj-2,knd)
          termd=dnewd*real(jj-2,knd)
110       m3temp=m3temp+term
          m3dtemp=m3dtemp+termd
          if((abs(term/m3temp)+abs(termd/m3dtemp)).lt.dcon) go to 130
          dold=dnew
          doldd=dnewd
120       continue
130     m3top=m3temp*cneun(l+1)
        im3e=ineue(l+1)
        iterm=log10(abs(m3top))
        im3e=im3e+iterm
        m3top=m3top*10.0e0_knd**(-iterm)
c        if(iopcs.eq.1) write(40,140) ndig,jneu,lim
c140     format(3x,'mc3/mc3d num. converged to ',i3,' digits in ',i10,
c     1         ' terms; ',i10,' available.')
c        if(iopcs.eq.2) write(40,150) ndig,jneu,lim
c150     format(3x,'ms3/ms3d num. converged to ',i3,' digits in ',i10,
c     1         ' terms; ',i10,' available.')
c
        m3c=m3top/m3bot
        if(iopcs.eq.2) m3c=m3c*con
        m3c=m3c*2.0e0_knd/pi
        if(ix.eq.0) m3c=-m3c
        if(iopcs.eq.1.and.ix.eq.1) m3c=-m3c
        iterm=log10(abs(m3c))
        im3e=im3e+iterm+iscale
        m3c=m3c*10.0e0_knd**(-iterm)
        if(abs(m3c).ge.1.0e0_knd) go to 160
        m3c=m3c*10.0e0_knd
        im3e=im3e-1
160	continue
        m3dtemp=m3dtemp*cneudr(l+1)
        if(iopcs.eq.2) m3dtemp=m3dtemp+m3temp/(xi*x2m1)
        m3dc=m3dtemp*cneun(l+1)/m3bot
        if(iopcs.eq.2) m3dc=m3dc*con
        m3dc=m3dc*2.0e0_knd/pi
        if(ix.eq.0) m3dc=-m3dc
        if(iopcs.eq.1.and.ix.eq.1) m3dc=-m3dc
        m3dc=m3dc*sqx2m1
        iterm=log10(abs(m3dc))
        im3de=iterm+ineue(l+1)+iscale
        m3dc=m3dc*10.0e0_knd**(-iterm)
        if(abs(m3dc).ge.1.0e0_knd) go to 170
        m3dc=10.0e0_knd*m3dc
        im3de=im3de-1
170     continue
        if(jneu.eq.lim) jneu=lim+20
        return
        end
c
c
         subroutine geteiga (l,cm,eiga1,eiga2,eiga3,eiga4,eiga5,eigaval)
c
c  purpose      : to calculate an estimate of the eigenvalue for
c                 c real. eigenvalues for c imaginary will be
c                 obtained from eigenvalues for c real using
c                 traditional relationships
c
c  parameters   :
c
c        input   : l        : l
c                  cm       : magnitude of c
c                eiga2-eiga5: previous eigenvalues
c
c        output  : eigaval  : estimate of the eigenvalue
c
        use param
c
c  real(knd) scalars
        real(knd) cm,d1,d2,d3,d4,eigaval,eiga1,eiga2,eiga3,eiga4,eiga5,
     1            phi,phi12,phi32,phi2,phi52,qu,qu2,qu3,qu4,qu6,w,w2,w3,
     2            w4,w5,w6,w7,w8,w9
c
        qu=cm*cm/4.0e0_knd
        qu2=qu*qu
        qu3=qu2*qu
        qu4=qu2*qu2
        qu6=qu4*qu2
c
c  after first four eigenvalues (l=0,1,2,3) have been computed, use
c  these via extrapolation to determine estimate for next eigenvalue
        if(l.gt.3) go to 30
c
c  use expansions in powers of qu for qu < 10, and
c  asymptotic expansions for qu >= 10.
        if(qu.ge.10.0e0_knd) go to 20
        if(qu.ge.1.0e0_knd) go to 10
c
c  estimates for qu < 1
        if(l.eq.0) eigaval=-0.5e0_knd*qu2+0.0546875e0_knd*qu4-
     1                      0.0125868e0_knd*qu6
        if(l.eq.1) eigaval=1.0e0_knd+qu-0.125e0_knd*qu2-0.0156e0_knd*qu3
        if(l.eq.2) eigaval=4.0e0_knd+0.417e0_knd*qu2-0.0552e0_knd*qu4
        if(l.eq.3) eigaval=9.0e0_knd+0.0625e0_knd*qu2+0.015625e0_knd*qu3
        go to 40
c  estimates for 1 < qu < 10
10      if(l.eq.0) eigaval=0.5543e0_knd-(0.883e0_knd*qu)-0.09639e0_knd*
     1                     qu2+0.004e0_knd*qu3
        if(l.eq.1) eigaval=0.81175e0_knd+1.3337e0_knd*qu-0.30892e0_knd*
     1                     qu2+0.01929e0_knd*qu3-0.0004946e0_knd*qu4
        if(l.eq.2) eigaval=3.3290504e0_knd+0.992e0_knd*qu-
     1                     0.0001829e0_knd*qu2-0.00866745e0_knd*qu3+
     2                     0.0003200972e0_knd*qu4
        if(l.eq.3) eigaval=8.9449274e0_knd-(0.1039356e0_knd*qu)+
     1                     0.190696e0_knd*qu2-0.01453021e0_knd*qu3+
     2                     0.000303573e0_knd*qu4
        go to 40
c  estimates for qu >= 10
20      w=l+l+1
        w2=w*w
        w3=w2*w
        w4=w2*w2
        w5=w3*w2
        w6=w3*w3
        w7=w3*w4
        w8=w4*w4
        w9=w5*w4
        phi=qu/w4
        phi12=sqrt(phi)
        phi2=phi*phi
        phi32=phi*phi12
        phi52=phi32*phi
        d1=5+34/w2+9/w4
        d2=33/w+410/w3+405/w5
        d3=63/w2+1260/w4+2943/w6+486/w8
        d4=527/w3+15617/w5+69001/w7+41607/w9
        eigaval=-2.0e0_knd*qu+2.0e0_knd*w*sqrt(qu)-((w2+1.0e0_knd)/
     1          8.0e0_knd)-((w+(3.0e0_knd/w))/(128.0e0_knd*phi12))-
     2          (d1/(4096.0e0_knd*phi))-(d2/(131072.0e0_knd*phi32))-
     3          (d3/(1048576.0e0_knd*phi2))-
     4          (d4/(33554432.0e0_knd*phi52))
        go to 40
c
c  third order extrapolation
30      eigaval=4.0e0_knd*eiga5-6.0e0_knd*eiga4+4.0e0_knd*eiga3-eiga2
        if(eigaval.lt.eiga5) eigaval=eiga5
40      eiga1=eiga2
        eiga2=eiga3
        eiga3=eiga4
        eiga4=eiga5
        return
        end
c
c
        subroutine convera (l,cm,limd,eiga1,eiga3,eiga4,eiga5,ndec,maxd,
     1                      enra,ienra,sgna,a01,ia01,blista,glista,
     2                      eigaval)
c
c  purpose      : to determine a converged eigenvalue using the
c                 boukwamp procedure. eigenvalues are obtained for
c                 c real; when c is imaginary, traditional relationships
c                 are used to obtain the necessary eigenvalues from
c                 those for c real.
c  parameters   :
c
c        input:   l      : l
c                 cm     : magnitude of c
c                 limd   : number of enra values computed
c                 eiga1  : eigenvalue for l-4
c                 eiga3  : eigenvalue for l-2
c                 eiga4  : eigenvalue for l-1
c                 ndec   : number of decimal digits available in real(knd)
c                          arithmetic
c                 maxd   : size of enr,blist,glist vectors
c                 eigaval: estimated value of the eigenvalue
c
c        output:  enra   : vector of coefficient ratios; enra(i) =
c                          a(sub 2i)/a(sub 2i-2) for l even or
c                          a(sub 2i+1)/a(sub 2i-1) for l odd
c                 ienra  : index i of last a coefficient ratio used
c                          in computing first term in the denominator
c                          of the eigenvalue correction
c                 sgna   : sign of the a coefficient a(sub l)
c                 a01    : characteristic of the ratio of the first a
c                          coefficent, either a(sub 0) if l is even or
c                          a(sub 1) if l is odd, to the a coefficient
c                          a(sub l)
c                 ia01   : exponent corresponding to a01
c                 eiga5  : converged eigenvalue for l
c
c        other:   blista : storage vector required by conver
c                 glista : storage vector required by conver
c
        use param
c
c  real(knd) scalars and vectors
        real(knd) a01,cm,cll,clu,cora,corb,de,dec,decinv,dl,
     1            dla,eiga1,eiga3,eiga4,eiga5,eigdec,eigaval,
     2            enririo,fl,qu,qu2,r,sgna
        real(knd) blista(maxd),glista(maxd),enra(maxd)
c
        qu=cm*cm/4.0e0_knd
        qu2=qu*qu
        dec=10.0e0_knd**(-ndec-1)
        decinv=10.0e0_knd**(ndec)
        eigdec=dec*100.0e0_knd
c
c  set the upper and lower bounds for the eigenvalue
        if(l.gt.0) cll=eiga4
        if(l.eq.1) clu=eigaval+(eigaval-eiga4)
        if(l.eq.2.or.l.eq.3) clu=eigaval+0.5e0_knd*(eigaval-eiga3)
        if(l.gt.3) clu=eigaval+0.5e0_knd*(eiga3-eiga1)
        l2=l/2
c
c  begin bouwkamp procedure
        fl=eigaval
        jnde=0
        ix=l-2*l2
        isc=2+ix
        j=1
c
c  compute the beta coeficients
          do 10 i=isc,limd,2
          blista(j)=qu2
          j=j+1
10        continue
        if(ix.eq.0) blista(1)=2.0e0_knd*qu2
        j=1
        id21=isc-1
        lim11=limd+1
        if(lim11.lt.id21) go to 30
c
c  compute the gamma coeficients
          do 20 i=id21,lim11,2
          r=real(i-1,knd)
          glista(j)=r*r
          j=j+1
20        continue
c  begin Bouwkamp procedure
30      ifc=1
        limc=limd/2-ix
        limdb=2*ienra+20
        if(l.eq.0) limdb=2*ienra
        if(limdb.gt.limd) limdb=limd
        lim2=limdb/2-ix 
        iglim=lim2+1
        irio=l2+1
        iw1=l2+2
40      enra(1)=eigaval-glista(1)
        if(ix.eq.1) enra(1)=enra(1)-qu
        if(l2.lt.1) go to 60
c
c  evaluate the continued fraction
          do 50 i=1,l2
          enra(i+1)=-blista(i)/enra(i)-glista(i+1)+eigaval
50        continue
60      enra(lim2)=-blista(lim2)/(glista(iglim)-eigaval)
        iw15=lim2-1
        ip=iw1+iw15
        if(iw15.lt.iw1) go to 80
          do 70 i=iw1,iw15
          ipi=ip-i
          enra(ipi)=-blista(ipi)/(glista(ipi+1)-eigaval+enra(ipi+1))
70        continue
80      enririo=-blista(irio)/(glista(irio+1)-eigaval+enra(irio+1))
        de=enririo*enririo/blista(irio)
        corb=de
        if(lim2.lt.iw1) go to 100
c
c  compute denominator
          do 90 i=iw1,lim2
          de=enra(i)*enra(i)/blista(i)*de
          corb=corb+de
          if(abs(de/corb).lt.dec) go to 100
90        continue
100     ienra=i
        if(ienra.lt.lim2-10) ienra=lim2-12
        cora=1.0e0_knd
        de=1.0e0_knd
        if(l2.lt.1) go to 120
          do 110 i=1,l2
          de=blista(irio-i)/(enra(irio-i)*enra(irio-i))*de
          cora=cora+de
          if(abs(de/cora).lt.dec) go to 120
110       continue
c
c  compute the correction to the eigenvalue
120     dl=(enririo-enra(irio))/(cora+corb)
        if(ifc.eq.1) dla=dl
        eigaval=dl+eigaval
c
c  eigenvalue accurate enough?
        if(abs(dl/eigaval).lt.eigdec) go to 130
        ifc=ifc+1
        if(ifc.lt.20) go to 40
130     continue
c
c  is the eigenvalue the correct one
c  if not then modify the original guess
        if(l.eq.0) go to 180
        if(eigaval.gt.cll) go to 140
c
c  converged to next lower eigenvalue of the same parity
        cll=fl
        go to 170
140     if(eigaval.lt.clu) go to 180
c
c  converged to the next higher eigenvalue of the same parity
        clu=fl
        go to 170
c
c  check to see if too many modifications are being made
c  if so then suspect error in the routine
170     jnde=jnde+1
        if(jnde.eq.50) go to 230
c
c  eigenvalue lies somewhere within the range established above
c  repeat procedure using the midpoint of this range
c  as the starting value
        eigaval=0.5e0_knd*(cll+clu)
        fl=eigaval
        ifc=1
        go to 40
c
180     continue
        eiga5=eigaval
c
c  compute enra using converged eigenvalue
          do i=2-ix,l2
          enra(i+1)=-blista(i)/enra(i)-glista(i+1)+eigaval
          end do
        enra(limc)=-blista(limc)/(glista(limc+1)-eigaval)
        iw15=limc-1
        ip=iw1+iw15
          if(iw15.ge.iw1) then
            do i=iw1,iw15
            ipi=ip-i
            enra(ipi)=-blista(ipi)/(glista(ipi+1)-eigaval+enra(ipi+1))
            end do
          end if
          if((l2.gt.2.and.abs(enra(l2+2)/enra(l2+3)).lt.0.01e0_knd).or.
     1        cm.lt.10.0e0_knd) then
          enra(l2+1)=-blista(l2+1)/(glista(l2+2)-eigaval+enra(l2+2))
          end if
        sgna=1.0e0_knd
        limit=limd/2-ix
          do 200 i=1,limit
          if(i.gt.l2) go to 190
          if(enra(i).lt.(0.0e0_knd))  sgna=-sgna
190       enra(i)=enra(i)/qu
200       continue
c
c  compute the ratio a(n=0)/a(n=l) for l even or a(n=1)/a(n=l) for l odd
c  this ratio has the characteristic a01 and the exponent ia01 
        ia01=0
   	a01=1.0e0_knd
        if(l2.eq.0) go to 220
          do 210 kjl=1,l2
          kkjl=l2-kjl+1
    	  a01=a01/enra(kkjl)
          if(abs(a01).gt.dec) go to 210
          a01=a01*decinv
          ia01=ia01-ndec
210       continue
        iterm=int(log10(abs(a01)))
        a01=a01*(10.0e0_knd**(-iterm))
        ia01=ia01+iterm
220     continue
        return
c
c  error printout
230     continue
c        write(20,240) l
c        write(30,240) l
        write(60,240) l
240     format(1x,'error in eigenvalue routine conver at l=',i5,/,
     1         1x,'this value may be inaccurate')
        eiga5=eigaval
        return
        end
c
c
        subroutine dnorma (l,limd,maxd,ndec,enra,ce0,
     1                     cepio2,cedpio2,jmfa,jsuba)
c
c  purpose    : to calculate several special values of the angular
c               functions ce and its first derivatives.
c
c  parameters :
c
c     input   : l       : l
c               limd    : equal to ~ twice the number of a coefficient
c                         ratios
c               maxd    : size of enra vector
c               ndec    : number of decimal digits available in
c                         real(knd) arithmetic
c               enra    : vector of a coefficient ratios
c
c     output  : ce0     : ce(0), divided by the a coefficient with
c                                subscript l
c               cepio2  : ce(pi/2), divided by the a coefficient with
c                                subscript l
c               cedpio2 : ce'(pi/2), divided by the a coefficient with
c                                subscript l
c               jmfa    : maximum number of a coefficient ratios
c                         required for convergence of all sums in
c                         this subroutine
c               jsuba   : subtraction error in calculating ce0
c
        use param
c
c  real(knd) scalars and vector
        real(knd) ce0,cedpio2,cepio2,dec,dnew,dold,sump,term
        real(knd) enra(maxd)
c
        dec=10.0e0_knd**(-ndec-2)
        l2=l/2
        ix=l-2*l2
        ix2=(l+1)/2-2*((l+1)/4)
        lim2=limd/2-ix
c
c  compute ce(0), divided by the a coefficient with subscript l
        ce0=1.0e0_knd
        sump=1.0e0_knd
        term=1.0e0_knd
        jlow=l+2
        jterm=l2
          do 30 j=jlow,limd,2
          jterm=jterm+1
          term=term*enra(jterm)
          if(term.gt.0.0e0_knd) sump=sump+term
          ce0=ce0+term
          if(abs(term).lt.dec) go to 40
30        continue
40      jlow=l
        jmfa=jterm
        if(jlow.lt.2) go to 60
        term=1.0e0_knd
        jterm=l2
          do 50 j=jlow,2,-2
          term=term/enra(jterm)
          if(term.gt.0.0e0_knd) sump=sump+term
          jterm=jterm-1
          ce0=ce0+term
          if(abs(term).lt.dec) go to 60
50        continue
60      continue
        jsuba=0
        if((sump*ce0).ne.0.0e0_knd) jsuba=
     1                             int(log10(abs(sump/ce0)))
        if(jsuba.lt.0) jsuba=0
        if(jsuba.gt.ndec) jsuba=ndec
        if(ce0.eq.0.0e0_knd.and.sump.ne.0.0e0_knd) jsuba=ndec
c        write(40,70) lim2,jmfa,jsuba
c70      format(3x,i10,' "a" coefficients; ce(0) converged in ',
c     1         i6,' terms; sub. error =',i3,' digits.')
        if(ix.eq.1) go to 125
c
c  compute ce(pi/2), divided by the a coefficient with subscript l
        term=1.0e0_knd
        cepio2=1.0e0_knd
          do 80 j=l2+1,lim2
          term=-term*enra(j)
          cepio2=cepio2+term
          if(abs(term/cepio2).lt.dec) go to 90
80        continue
90      continue
        jterm=j
        jmfa=max(jmfa,jterm)
        if(l2.lt.1) go to 110
        term=1.0e0_knd
          do 100 j=l2,1,-1
          term=-term/enra(j)
          cepio2=cepio2+term
          if(abs(term/cepio2).lt.dec) go to 110
100       continue
110     nsub=0
        if(ix2.eq.1) cepio2=-cepio2
125     continue
        if(ix.eq.0) go to 180
c
c  compute ce'(pi/2), divided by the a coefficient with subscript l
        dold=1.0e0_knd
        cedpio2=real(l,knd)
          do 130 j=l2+1,lim2
          dnew=-dold*enra(j)
          term=dnew*real(j+j+ix,knd)
          cedpio2=cedpio2+term
          if(abs(term/cedpio2).lt.dec) go to 140
          dold=dnew
130       continue
140     continue
        jterm=j
        jmfa=max(jmfa,jterm)
        if(l2.lt.1) go to 160
        dold=1.0e0_knd
          do 150 j=l2,1,-1
          dnew=-dold/enra(j)
          term=dnew*real(j+j+ix-2,knd)
          cedpio2=cedpio2+term
          if(abs(term/cedpio2).lt.dec) go to 160
          dold=dnew
150       continue
160     if(ix2.eq.1) cedpio2=-cedpio2
180     continue
        return
        end
c
c
        subroutine geteigb (l,cm,eigb1,eigb2,eigb3,eigb4,eigb5,eigbval)
c
c  purpose      : to calculate an estimate of the eigenvalue for
c                 c real. eigenvalues for c imaginary will be
c                 obtained from eigenvalues for c real using
c                 traditional relationships
c
c  parameters   :
c
c        input   : l        : l
c                  cm       : cm
c                eigb2-eigb5: previous eigenvalues
c
c        output  : eigbval  : estimate of the eigenvalue
c
        use param
c
c  real(knd) scalars
        real(knd) cm,d1,d2,d3,d4,eigbval,eigb1,eigb2,eigb3,eigb4,eigb5,
     1            phi,phi12,phi32,phi2,phi52,qu,qu2,qu3,qu4,qu6,w,w2,w3,
     2            w4,w5,w6,w7,w8,w9
c
        qu=cm*cm/4.0e0_knd
        qu2=qu*qu
        qu3=qu2*qu
        qu4=qu2*qu2
        qu6=qu4*qu2
c
c  after first four eigenvalues (l=0,1,2,3) have been computed, use
c  these via extrapolation to determine estimate for next eigenvalue
        if(l.gt.4) go to 30
c
c  use expansions in powers of qu for qu < 10, and
c  asymptotic expansions for qu >= 10.
        if(qu.ge.10.0e0_knd) go to 20
        if(qu.ge.1.0e0_knd) go to 10
c
c  estimates for qu < 1
        if(l.eq.1) eigbval=1.0e0_knd-qu-0.125e0_knd*qu2+0.0156e0_knd*qu3
        if(l.eq.2) eigbval=4.0e0_knd+0.08333333e0_knd*qu2+
     1                     (3.6169e-4_knd)*qu4
        if(l.eq.3) eigbval=9.0e0_knd+0.0625e0_knd*qu2-0.015625e0_knd*qu3
        if(l.eq.4) eigbval=16.0e0_knd+0.033333333e0_knd*qu2-
     1                     (3.66898e-04_knd)*qu4
        go to 40
c  estimates for 1 < qu < 10
10      if(l.eq.1) eigbval=1.10427e0_knd-(1.152218e0_knd*qu)-
     1                     0.0548246e0_knd*qu2+
     2                     0.0019711e0_knd*qu3
        if(l.eq.2) eigbval=4.00909e0_knd-(0.004732542e0_knd*qu)-
     1                     0.08725329e0_knd*qu2
     2                     +0.00238446e0_knd*qu3
        if(l.eq.3) eigbval=8.771735e0_knd+0.268987e0_knd*qu-
     1                     0.03569325e0_knd*qu2+
     2                     0.00009369364e0_knd*qu3
        if(l.eq.4) eigbval=15.744e0_knd+0.1907483e0_knd*qu+
     1                     0.0038216144e0_knd*qu2-
     2                     0.000708719e0_knd*qu3
        go to 40
c  estimates for qu >= 10
20      w=l+l-1
        w2=w*w
        w3=w2*w
        w4=w2*w2
        w5=w3*w2
        w6=w3*w3
        w7=w3*w4
        w8=w4*w4
        w9=w5*w4
        phi=qu/w4
        phi12=sqrt(phi)
        phi2=phi*phi
        phi32=phi*phi12
        phi52=phi32*phi
        d1=5+34/w2+9/w4
        d2=33/w+410/w3+405/w5
        d3=63/w2+1260/w4+2943/w6+486/w8
        d4=527/w3+15617/w5+69001/w7+41607/w9
        eigbval=-2.0e0_knd*qu+2.0e0_knd*w*sqrt(qu)-((w2+1.0e0_knd)/
     1          8.0e0_knd)-((w+(3.0e0_knd/w))/(128.0e0_knd*phi12))-
     2          (d1/(4096.0e0_knd*phi))-(d2/(131072.0e0_knd*phi32))-
     3          (d3/(1048576.0e0_knd*phi2))-(d4/(33554432.0e0_knd*
     4          phi52))
        go to 40
c
c  third order extrapolation
30      eigbval=4.0e0_knd*eigb5-6.0e0_knd*eigb4+4.0e0_knd*eigb3-eigb2
        if(eigbval.lt.eigb5) eigbval=eigb5
40      eigb1=eigb2
        eigb2=eigb3
        eigb3=eigb4
        eigb4=eigb5
        return
        end
c
c
        subroutine converb (l,cm,limd,eigb1,eigb3,eigb4,eigb5,ndec,maxd,
     1                      enrb,ienrb,sgnb,b12,ib12,blistb,glistb,
     2                      eigbval)
c
c  purpose      : to determine a converged eigenvalue using the
c                 boukwamp procedure for c real. eigenvalues for
c                 c imaginary will be obtained from eigenvalues
c                 for c real using traditional relationships
c  parameters   :
c
c        input:   l      : l
c                 cm     : magnitude of c
c                 eigb1  : eigenvalue for l-4
c                 eigb3  : eigenvalue for l-2
c                 eigb4  : eigenvalue for l-1
c                 ndec   : number of decimal digits available in
c                          real(knd) arithmetic
c                 maxd   : size of enr,blist,glist vectors
c                 eigbval: estimated value of the eigenvalue
c
c        output:  enrb   : vector of coefficient ratios; enrb(i) =
c                          b(sub 2i)/b(sub 2i-2) for l even or
c                          b(sub 2i+1)/b(sub 2i-1) for l odd; Since
c                          b(sub 0) is not used, enrb(1) is not defined
c                          when l is even 
c                 ienrb  : index i of last b coefficient ratio used
c                          in computing first term in the denominator
c                          of the eigenvalue correction
c                 sgnb   : sign of the b coefficient b(sub l)
c                 b01    : characteristic of the ratio of the first b
c                          coefficent, either b(sub 2) if l is even or
c                          b(sub 1) if l is odd, to the b coefficient
c                          b(sub l)
c
c        other:   blistb : storage vector required by conver
c                 glistb : storage vector required by conver
c
        use param
c
c  real(knd) scalars and vectors
        real(knd) b12,cm,cll,clu,cora,corb,de,dec,decinv,dl,
     1            dla,eigb1,eigb3,eigb4,eigb5,eigdec,eigbval,
     2            enrc,fl,qu,qu2,r,sgnb
        real(knd) blistb(maxd),glistb(maxd),enrb(maxd)
c
        qu=cm*cm/4.0e0_knd
        qu2=qu*qu
        l2=l/2
        ix=l-2*l2
        dec=10.0e0_knd**(-ndec-1)
        decinv=10.0e0_knd**(ndec)
        eigdec=dec*100.0e0_knd
c
c  set the upper and lower bounds for the eigenvalue
        if(l.gt.1) cll=eigb4
        if(l.eq.2) clu=eigbval+(eigbval-eigb4)
        if(l.eq.3.or.l.eq.4) clu=eigbval+0.5e0_knd*(eigbval-eigb3)
        if(l.gt.4) clu=eigbval+0.5e0_knd*(eigb3-eigb1)
c
c  begin bouwkamp procedure
        fl=eigbval
        jnde=0
        isc=2+ix
        j=1
c
c  compute the beta coeficients
          do 10 i=isc,limd,2
          blistb(j)=qu2
          j=j+1
10        continue
        j=1
        id21=isc-1
        lim11=limd+1
        if(lim11.lt.id21) go to 30
c
c  compute the gamma coeficients
          do 20 i=id21,lim11,2
          r=real(i-1,knd)
          glistb(j)=r*r
          j=j+1
20        continue
30      ifc=1
        limc=limd/2-ix
        limdb=2*ienrb+20
        if(l.eq.1) limdb=2*ienrb
        if(limdb.gt.limd) limdb=limd
        lim2=limdb/2-ix
        iglim=lim2+1
        irio=l2+1
        iw1=l2+2
40      if(ix.eq.1) enrb(1)=eigbval-glistb(1)+qu
        if(ix.eq.0) enrb(2)=eigbval-glistb(2)
        if(l2.lt.(2-ix)) go to 60
c
c  evaluate the continued fraction
          do 50 i=2-ix,l2
          enrb(i+1)=-blistb(i)/enrb(i)-glistb(i+1)+eigbval
50        continue
60      enrb(lim2)=-blistb(lim2)/(glistb(iglim)-eigbval)
        iw15=lim2-1
        ip=iw1+iw15
        if(iw15.lt.iw1) go to 80
          do 70 i=iw1,iw15
          ipi=ip-i
          enrb(ipi)=-blistb(ipi)/(glistb(ipi+1)-eigbval+enrb(ipi+1))
70        continue
80      enrc=-blistb(irio)/(glistb(irio+1)-eigbval+enrb(irio+1))
        de=enrc*enrc/blistb(irio)
        corb=de
        if(lim2.lt.iw1) go to 100
c
c  compute denominator
          do 90 i=iw1,lim2
          de=enrb(i)*enrb(i)/blistb(i)*de
          corb=corb+de
          if(abs(de/corb).lt.dec) go to 100
90        continue
100     ienrb=i
        if(ienrb.lt.lim2-10) ienrb=lim2-12
        cora=1.0e0_knd
        de=1.0e0_knd
        if(l2+ix.lt.2) go to 120
          do 110 i=1,l2-1+ix
          de=blistb(irio-i)/(enrb(irio-i)*enrb(irio-i))*de
          cora=cora+de
          if(abs(de/cora).lt.dec) go to 120
110       continue
c
c  compute the correction to the eigenvalue
120     dl=(enrc-enrb(irio))/(cora+corb)
        if(ifc.eq.1) dla=dl
        eigbval=dl+eigbval
c
c  eigenvalue accurate enough?
        if(abs(dl/eigbval).lt.eigdec) go to 130
        ifc=ifc+1
        if(ifc.lt.20) go to 40
130     continue
c
c  is the eigenvalue the correct one
c  if not then modify the original guess
        if(l.lt.2) go to 180
        if(eigbval.gt.cll) go to 140
c
c  converged to next lower eigenvalue of the same parity
        cll=fl
        go to 170
140     if(eigbval.lt.clu) go to 180
c
c  converged to the next higher eigenvalue of the same parity
        clu=fl
        go to 170
c
c  check to see if too many modifications are being made
c  if so then suspect error in the routine
170     jnde=jnde+1
        if(jnde.eq.50) go to 230
c
c  eigenvalue lies somewhere within the range established above
c  repeat bouwkamp procedure using the midpoint value as the
c  starting value
        eigbval=.5e0_knd*(cll+clu)
        fl=eigbval
        ifc=1
        go to 40
c
180     continue
        eigb5=eigbval
          do i=2-ix,l2
          enrb(i+1)=-blistb(i)/enrb(i)-glistb(i+1)+eigbval
          end do
        enrb(limc)=-blistb(limc)/(glistb(limc+1)-eigbval)
        iw15=limc-1
        ip=iw1+iw15
          if(iw15.ge.iw1) then
            do i=iw1,iw15
            ipi=ip-i
            enrb(ipi)=-blistb(ipi)/(glistb(ipi+1)-eigbval+enrb(ipi+1))
            end do
          end if
          if((l2.gt.2.and.abs(enrb(l2+2)/enrb(l2+3)).lt.0.01e0_knd).or.
     1        cm.lt.10.0e0_knd) then
          enrb(l2+1)=-blistb(l2+1)/(glistb(l2+2)-eigbval+enrb(l2+2))
          end if
        sgnb=1.0e0_knd
        jlim=1
        if(ix.eq.0) jlim=2
          do 200 j=jlim,limc
          if(j.gt.l2) go to 190
          if(enrb(j).lt.(0.0e0_knd))  sgnb=-sgnb
190       enrb(j)=enrb(j)/qu
200       continue
        if(ix.eq.0) enrb(1)=0.0e0_knd
c
c  compute b1(cm,l) or b2(cm,l)
        ib12=0
   	b12=1.0e0_knd
        if(l2.eq.0) go to 220
        kjlupp=l2
        if(ix.eq.0) kjlupp=l2-1
        if(kjlupp.eq.0) go to 220
          do 210 kjl=1,kjlupp
          kkjl=l2-kjl+1
    	  b12=b12/enrb(kkjl)
          if(abs(b12).gt.dec) go to 210
          b12=b12*decinv
          ib12=ib12-ndec
210       continue
        iterm=int(log10(abs(b12)))
        b12=b12*(10.0e0_knd**(-iterm))
        ib12=ib12+iterm
220     continue
        return
c
c  error printout
230     continue
c        write(20,240) l
c        write(30,240) l
        write(60,240) l 
240     format(1x,'error in eigenvalue routine conver at l=',i5,/,
     1         1x,'this value may be inaccurate')
        eigb5=eigbval
        return
        end
c
c
        subroutine dnormb (l,limd,maxd,ndec,enrb,sed0,sepio2,
     1                     sedpio2,jmfb,jsubb)
c
c  purpose    : to calculate several special values of the angular
c               functions se and its first derivative.
c
c  parameters :
c
c     input   : l       : l
c               limd    : equal to ~twice the number of b coefficient
c                         ratios
c               maxd    : size of enra vector
c               ndec    : number of decimal digits available in
c                         real(knd) arithmetic
c               enrb    : vector of b coefficient ratios
c
c     output:   sed0    : se'(0), divided by the b coefficient with
c                                subscript l
c               sepio2  : se(pi/2), divided by the b coefficient with
c                                subscript l
c               sedpio2 : se'(pi/2), divided by the b coefficient with
c                                subscript l
c               jmfb    : maximum number of b coefficient ratios
c                         required for convergence of all sums in
c                         this subroutine
c               jsubb   : subtraction error in calculating sed0
c
        use param
c
c  real(knd) scalars and vector
        real(knd) dec,dnew,dold,sed0,sedpio2,sepio2,sump,term
        real(knd) enrb(maxd)
c
        dec=10.0e0_knd**(-ndec-2)
        l2=l/2
        ix=l-2*l2
        ix2=(l+1)/2-2*((l+1)/4)
        lim2=limd/2-ix
        jlim=1
        if(ix.eq.0) jlim=2
c
c  compute se'(0), divided by the b coefficient with subscript l
        sed0=real(l,knd)
        sump=sed0
        dold=1.0e0_knd
          do 30 j=l2+1,lim2
          dnew=dold*enrb(j)
          term=dnew*real(j+j+ix,knd)
          if(term.gt.0.0e0_knd) sump=sump+term
          sed0=sed0+term
          if(abs(term/sed0).lt.dec) go to 40
          dold=dnew
30        continue
40      jmfb=j
        if(l2.lt.jlim) go to 60
        dold=1.0e0_knd
          do 50 j=l2,jlim,-1
          dnew=dold/enrb(j)
          term=dnew*real(j+j+ix-2,knd)
          if(term.gt.0.0e0_knd) sump=sump+term
          sed0=sed0+term
          if(abs(term/sed0).lt.dec) go to 60
          dold=dnew
50        continue
60      continue
        jsubb=0
        if((sump*sed0).ne.0.0e0_knd) jsubb=
     1                               int(log10(abs(sump/sed0)))
        if(jsubb.lt.0) jsubb=0
        if(jsubb.gt.ndec) jsubb=ndec
        if(sed0.eq.0.0e0_knd.and.sump.ne.0.0e0_knd) jsubb=ndec
c        write(40,70) lim2,jmfb,jsubb
c70      format(3x,i10,' "b" coefficients; sed(0) converged in',
c     1         i6,' terms; sub. error =',i3,' digits.')
        if(ix.eq.0) go to 125
c
c  compute se(pi/2), divided by the b coefficient with subscript l
        term=1.0e0_knd
        sepio2=1.0e0_knd
          do 80 j=l2+1,lim2
          term=-term*enrb(j)
          sepio2=sepio2+term
          if(abs(term/sepio2).lt.dec) go to 90
80        continue
90      continue
        jterm=j
        jmfb=max(jmfb,jterm)
        if(l2.lt.1) go to 110
        term=1.0e0_knd
          do 100 j=l2,1,-1
          term=-term/enrb(j)
          sepio2=sepio2+term
          if(abs(term/sepio2).lt.dec) go to 110
100       continue
110     continue
        if(ix2.eq.0) sepio2=-sepio2
125     continue
        if(ix.eq.1) go to 180
c
c  compute se'(pi/2), divided by the b coefficient with subscript l
        dold=1.0e0_knd
        sedpio2=real(l,knd)
          do 130 j=l2+1,lim2
          dnew=-dold*enrb(j)
          term=dnew*real(j+j+ix,knd)
          sedpio2=sedpio2+term
          if(abs(term/sedpio2).lt.dec) go to 140
          dold=dnew
130       continue
140     continue
        jterm=j
        jmfb=max(jmfb,jterm)
        if(l2.lt.2) go to 160
        dold=1.0e0_knd
          do 150 j=l2,2,-1
          dnew=-dold/enrb(j)
          term=dnew*real(j+j+ix-2,knd)
          sedpio2=sedpio2+term
          if(abs(term/sedpio2).lt.dec) go to 160
          dold=dnew
150       continue
160     if(ix2.eq.1) sedpio2=-sedpio2
180     continue
        return
        end
c
c
        subroutine sincos(limsc,maxp,narg,isq,barg,sine,cosi)
c
        use param
c
c  real(knd) scalars and vectors
        real(knd) cosi1,sine1
        real(knd) barg(narg),sine(narg,maxp),cosi(narg,maxp)
c
          do 40 i=1,narg
          sine(i,1)=0.0e0_knd
          cosi(i,1)=1.0e0_knd
          if(isq.eq.-1) go to 10
          sine1=sin(barg(i))
          cosi1=cos(barg(i))
          go to 20
10        cosi1=sin(barg(i))
          sine1=cos(barg(i))
20          do 30 j=2,limsc
            sine(i,j)=cosi1*sine(i,j-1)+sine1*cosi(i,j-1)
            cosi(i,j)=cosi1*cosi(i,j-1)-sine1*sine(i,j-1)
30          continue
40        continue
        return
        end
c
c
        subroutine cbesj (cm,x,limj,maxj,numlp,maxlp,limjt,ndec,limbes,
     1                    cbesf,cbesn,ibese,cbesdf,cbesdr,iopr)
c
c  purpose    : to calculate ratios of cylindrical Bessel functions of
c               the first kind and ratios of their first derivatives
c               for given c and xi. to calculate the Bessel function
c               characteristics and exponents and the ratios of the
c               first derivatives to the functions.
c
c  parameters :
c
c     input     cm     : c, since c is real here
c               x      : value of radial argument equal to xi or
c                        sqrt(xi*xi-1) or 0.5/[xi + sqrt(xi*xi-1)]
c                        depending on the expansion used to
c                        calculate the radial functions. the
c                        argument of the Bessel functions is cm*x.
c               limj   : the number of cylindrical Bessel function
c                        ratios calculated for given lnum, c, and ndec.
c               maxj   : size of cbesf and sbesdf vectors
c               numlp  : the number of values of scale factors
c                        and ratios of first derivatives to functions
c                        that are calculated
c               maxlp  : dimension of vectors of scale factors and
c                        ratios of first derivatives to functions
c               limjt  : Bessel function order to begin backward
c                        recursion
c               ndec   : number of decimal digits for real(knd)
c                        arithmetic
c               limbes : 2*c*xi+2*ndec
c
c     output  : cbesf  : ratios of cylindrical Bessel functions of the
c                        first kind
c               cbesn  : characteristics for the cylindrical
c                        Bessel functions
c               ibese  : exponents for the cylindrical
c                        Bessel functions
c               cbesdf : ratios of first derivatives of the cylindrical
c                        Bessel functions
c               cbesdr : ratios of first derivatives of the cylindrical
c                        Bessel functions to the corresponding Bessel
c                        functions
c
        use param
c
c  real(knd) scalars and vectors
        real(knd) besj0,cm,cx,dec,rn,sumj,temp0,temp1,term,x
        real(knd) cbesdf(maxj),cbesdr(maxj),cbesf(maxj),cbesn(maxlp),
     1            ratio(limjt)
c
c  integer vector
        dimension ibese(maxlp)
c
        cx=cm*x
        dec=10.0e0_knd**(-ndec-2)
        lim1=max(limj+ndec,limbes)
c
c  compute first kind Bessel function ratios
c        cbesf(k)= J(n=k,c*x)/J(n=k-1,c*x)
c        cbesn(k)= (J(n=k-1),c*x))*10.0e0_knd**(-ibese(k))
c
c  use backward recursion to compute                                                                                                                                 ratios:
c       J(n,c*x)/J(n-1,c*x) = 1/[(2*n)/(c*x) - J(n+1,c*x)/J(n,c*x)]
c
        temp0=cx/(2.0e0_knd*real(lim1+1,knd))
        if(limjt.ge.lim1) go to 15
          do 10 n=lim1,limjt,-1
          rn=real(n+n,knd)
          temp1=1.0e0_knd/(rn/cx-temp0)
          temp0=temp1
10        continue
15      if(limjt.ge.lim1) temp1=cx/(2.0e0_knd*real(limjt,knd))
        ratio(limjt)=temp1
          do 20 n=limjt-1,1,-1
          rn=real(n+n,knd)
          ratio(n)=1.0e0_knd/(rn/cx-ratio(n+1))
20        continue
          do 30 n=1,limj
          cbesf(n)=ratio(n)
30        continue
c
c  calculate the zeroth order function J(0,c*x) by use of the ratios
c  calculated above and the formula
c
c       1 = J(0,c*x) + 2J(2,c*x) + 2J(4,c*x) + ..., or
c
c       J(0,c*x) = 1/[1 + 2(J(2,c*x)/J(0,c*x)) + 2(J(4,c*x)/J(0,c*x)) +
c                  ...]
c
        term=2.0e0_knd
        sumj=1.0e0_knd
          do 40 j=2,limbes,2
          term=term*ratio(j-1)*ratio(j)
          sumj=sumj+term
          if(abs(term/sumj).lt.dec) go to 50
40        continue
50      continue
        besj0=1.0e0_knd/sumj
c
c  calculate the characteristics and exponents for the
c  Bessel functions by forward operation on the function
c  ratios starting with the calculated value for J(0,c*x).
        ibese(1)=int(log10(abs(besj0)))
        cbesn(1)=besj0*10.0e0_knd**(-ibese(1))
        cbesn(2)=cbesn(1)*cbesf(1)
        ibese(2)=int(log10(abs(cbesn(2))))
        cbesn(2)=cbesn(2)*10.0e0_knd**(-ibese(2))
        ibese(2)=ibese(2)+ibese(1)
          do 60 n=3,numlp
          cbesn(n)=cbesn(n-1)*cbesf(n-1)
          ibese(n)=int(log10(abs(cbesn(n))))
          cbesn(n)=cbesn(n)*10.0e0_knd**(-ibese(n))
          ibese(n)=ibese(n)+ibese(n-1)
60        continue
c
c  calculate the ratios of the first derivatives of successive
c  Bessel functions using corresponding function ratios
c  calculate the ratios of the first derivative to the corresponding
c  cylindrical Bessel function
          do 70 n=1,limj
          rn=real(n-1,knd)
          cbesdf(n)=(cx-(rn+1.0e0_knd)*cbesf(n))/(rn-cx*cbesf(n))
          cbesdr(n)=(rn/cx)-cbesf(n)
70        continue
        if(iopr.eq.1) go to 90
c
c  convert ratios of successive functions or derivatives to ratios of
c  functions or derivatives of the same parity
          do 80 n=limj,2,-1
          cbesf(n)=cbesf(n-1)*cbesf(n)
80        cbesdf(n)=cbesdf(n-1)*cbesdf(n)
90      continue
        return
        end
c
c
        subroutine cbesi (cm,x,limj,maxj,numlp,maxlp,limjt,ndec,limbes,
     1                    cbesf,cbesn,ibese,cbesdf,cbesdr,iopr,nex)
c
c  purpose    : to calculate ratios of modified cylindrical Bessel
c               functions of the first kind and ratios of their first
c               derivatives for given c and x. to calculate the
c               modified Bessel function characteristics and exponents
c               and ratios of the first derivatives to the functions.
c
c  parameters :
c
c     input     cm     : mgnitude of c, c being positive imaginary
c               x      : value of radial argument for modified Bessel
c                        functions equal to cm*xi or cm*sqrt(xi*xi-1) or
c                        cm*0.5*[xi - sqrt(xi*xi-1)] or cm*0.5*[xi +
c                        sqrt(xi*xi-1)] depending on the expansion
c                        used to calculate the radial functions
c               limj   : the number of cylindrical Bessel function
c                        ratios calculated for given lnum, c, and ndec.
c               maxj   : size of cbesf and cbesdf vectors
c               numlp  : the number of values of scale factors
c                        and ratios of first derivatives to functions
c                        that are calculated
c               maxlp  : dimension of vectors of scale factors and
c                        ratios of first derivatives to functions
c               limjt  : Bessel function order to begin backward
c                        recursion
c               ndec   : number of decimal digits for real(knd)
c                        arithmetic
c               limbes : 2*x+2*ndec
c
c     output  : cbesf  : ratios of modified cylindrical Bessel functions
c                        of the first kind
c               cbesn  : characteristics for the modified cylindrical
c                        Bessel functions
c               ibese  : exponents for the modified cylindrical
c                        Bessel functions
c               cbesdf : ratios of derivatives of the modified
c                        cylindrical Bessel functions
c               cbesdr : ratios of the derivatives of the modified
c                        cylindrical Bessel functions to their
c                        corresponding Bessel functions
c
        use param
c
c  real(knd) scalars and vectors
        real(knd) besi0,cm,cx,dec,fac,rn,sumj,temp0,temp1,term,x
        real(knd) cbesdf(maxj),cbesdr(maxj),cbesf(maxj),cbesn(maxlp),
     1            ratio(limjt)
c
c  integer vector
        dimension ibese(maxlp)
c
        cx=cm*x
        dec=10.0e0_knd**(-ndec-2)
        lim1=max(limj+ndec,limbes)
c
c  compute first kind modified Bessel function ratios
c        cbesf(k)= I(n=k,cm*x)/I(n=k-1,cm*x)
c        cbesn(k)= (I(n=k-1),cm*x))*10.0e0_knd**(-ibese(k))
c
c  use backward recursion to compute                                                                                                                                 ratios:
c       I(n,cm*x)/I(n-1,cm*x) = 1/((2*n)/(cm*x) + I(n+1,cm*x)/I(n,cm*x))
c
        temp0=cx/(2.0e0_knd*real(lim1+1,knd))
        if(limjt.ge.lim1) go to 20
          do 10 n=lim1,limjt,-1
          rn=real(n+n,knd)
          temp1=1.0e0_knd/(rn/cx+temp0)
          temp0=temp1
10        continue
20      if(limjt.ge.lim1) temp1=cx/(2.0e0_knd*real(limjt,knd))
        ratio(limjt)=temp1
          do 30 n=limjt-1,1,-1
          rn=real(n+n,knd)
          ratio(n)=1.0e0_knd/(rn/cx+ratio(n+1))
30        continue
          do 40 n=1,limj
          cbesf(n)=ratio(n)
40        continue
c
c  calculate the zeroth order function I(0,cm*x) divided by
c  exp(cm*x) by use of the ratios calculated above and the formula
c
c       exp(cm*x) = I(0,cm*x) + 2I(1,cm*x) + 2I(2,cm*x) + ..., or
c
c       I(0,cm*x) = exp(cm*x)/[1 + 2(I(1,cm*x)/I(0,cm*x)) +
c                  2(I(2,cm*x)/I(0,cm*x)) + ...]
c
        term=2.0e0_knd
        sumj=1.0e0_knd
          do 50 j=1,limbes
          term=term*ratio(j)
          sumj=sumj+term
          if(abs(term/sumj).lt.dec) go to 60
50        continue
60      continue
        besi0=1.0e0_knd/sumj
c
c  for all cm*x, calculate the amplitude and sign scale
c  factors by forward operation on the Bessel function
c  ratios starting with the calculated value for i(0,cm*x)
c  divided by exp(cm*x).
        ncx=2*int(cx*log10(exp(1.0e0_knd)))/nex
          if(ncx.lt.1) then
          cbesn(1)=besi0*exp(cx)
          iterm=int(log10(abs(cbesn(1))))
          cbesn(1)=cbesn(1)*(10.0e0_knd**(-iterm))
          ibese(1)=iterm
          else
          fac=exp(cx/real(ncx,knd))
          cbesn(1)=besi0*fac
          ibese(1)=int(log10(abs(cbesn(1))))
          cbesn(1)=cbesn(1)*(10.0e0_knd**(-ibese(1)))
            do j=2,ncx
            cbesn(1)=cbesn(1)*fac
            iterm=int(log10(abs(cbesn(1))))
            cbesn(1)=cbesn(1)*(10.0e0_knd**(-iterm))
            ibese(1)=ibese(1)+iterm
            end do
          end if
        cbesn(2)=cbesn(1)*cbesf(1)
        ibese(2)=int(log10(abs(cbesn(2))))
        cbesn(2)=cbesn(2)*10.0e0_knd**(-ibese(2))
        ibese(2)=ibese(2)+ibese(1)
          do 70 n=3,numlp
          cbesn(n)=cbesn(n-1)*cbesf(n-1)
          ibese(n)=int(log10(abs(cbesn(n))))
          cbesn(n)=cbesn(n)*10.0e0_knd**(-ibese(n))
          ibese(n)=ibese(n)+ibese(n-1)
70        continue
c
c  calculate the ratios of the first derivatives of successive
c  modified Bessel functions using corresponding function ratios
c  calculate the ratios of the first derivative to the corresponding
c  modified cylindrical Neumann function

          do 80 n=1,limj
          rn=real(n-1,knd)
          cbesdf(n)=(cx-(rn+1.0e0_knd)*cbesf(n))/(rn+cx*cbesf(n))
          cbesdr(n)=(rn/cx)+cbesf(n)
80        continue
c
        if(iopr.eq.1) go to 100
c
c  convert ratios of successive functions and derivatives to ratios of
c  functions and derivatives of the same parity
          do 90 n=limj,2,-1
          cbesf(n)=cbesf(n-1)*cbesf(n)
          cbesdf(n)=cbesdf(n-1)*cbesdf(n)
90        continue
100     continue
        return
        end
c
c
        subroutine cneuy (cm,x,limn,maxn,numlp,maxlp,ndec,limbes,pi,
     1                    gamma,cneuf,cneun,ineue,cneudf,cneudr,iopr)
c
c  purpose    : to calculate same parity ratios of cylindrical Bessel
c               functions of the second kind (Neumann functions) and
c               same parity ratios of their first derivatives for
c               given c and x.
c               to calculate the Neumann function and Neumann function
c               characteristics and exponents.
c               to calculate the ratios of the first derivatives to
c               corresponding cylindrical Neumann functions
c
c  parameters :
c
c     input     cm     : c, since c is real here
c               x      : x
c               limn   : the number of cylindrical Neumann function
c                        ratios calculated for given lnum, c, and ndec
c               maxn   : size of cneuf, cneudf and cneudr vectors
c               numlp  : the number of values of scale factors
c                        that are calculated
c               maxlp  : dimension of vectors of scale factors
c               ndec   : number of decimal digits for real(knd)
c                        arithmetic
c               limbes : 2*c*x+2*ndec
c
c     output  : cneuf  : same parity ratios of cylindrical Neumann
c                        functions
c               cneun  : characteristic for the cylindrical
c                        Neumann functions
c               ineue  : exponent for the cylindrical
c                        Neumann functions
c               cneudf : same parity ratios of first derivatives
c                        of cylindrical Neumann functions
c               cneudr : ratios of first derivatives to the
c                        corresponding cylindrical Neumann functions
c
        use param
c
c  real(knd) scalars and vectors
        real(knd) besj0,cm,cx,cy0,cy1,dec,gamma,pi,rn,sumj,term,x
        real(knd) cbesf(limbes),cneuf(maxn),cneun(maxlp),cneudf(maxn),
     1            cneudr(maxn)
c
c  integer vectors
        dimension ineue(maxlp)
        cx=cm*x
        dec=10.0e0_knd**(-ndec-2)
c  calculate ratios of successive cylindrical Bessel functions of the
c  first kind
        cbesf(limbes)=cx/(2.e0_knd*real(limbes,knd))
          do 10 n=limbes-1,1,-1
          rn=real(n+n,knd)
          cbesf(n)=1.0e0_knd/(rn/cx-cbesf(n+1))
10        continue
c
c  calculate the zeroth order function J(n=0,c*x) by use of the ratios
c  calculated above and the formula
c
c       1 = J(0,c*x) + 2J(2,c*x) + 2J(4,c*x) + ..., or
c
c       J(0,c*x) = 1/[1 + 2(J(2,c*x)/J(0,c*x)) + 2(J(4,c*x)/J(0,c*x)) +
c                  ...]
c
        term=2.0e0_knd
        sumj=1.0e0_knd
          do 20 j=2,limbes,2
          term=term*cbesf(j-1)*cbesf(j)
          sumj=sumj+term
          if(abs(term/sumj).lt.dec) go to 30
20        continue
30      continue
        besj0=1.0e0_knd/sumj
c
c  calculate cy0 = Y(n=0,c*x) using J(n=0,c*x) and J ratios
c
        term=besj0*cbesf(1)*cbesf(2)
        sumj=term
          do 40 j=2,limbes/2
          term=-term*cbesf(j+j-1)*cbesf(j+j)
          sumj=sumj+term/real(j,knd)
          if(abs(term/sumj).lt.dec) go to 50
40        continue
50      cy0=2.0e0_knd*((log(cx/2.0e0_knd)+gamma)*besj0+2.0e0_knd*sumj)
     1            /pi
        cy1=(cy0*besj0*cbesf(1)-2.0e0_knd/(pi*cx))/besj0
c
c  compute ratios of successive cylindrical Neumann functions
c  using forward recursion from cy1/cy0
c
c       y(n+1,c*x)/y(n,c*x)=(2*n)/(c*x)-1/(y(n,c*x)/y(n-1,c*x))
c
        cneuf(1)=cy1/cy0
          do 60 n=1,limn-1
          rn=real(n+n,knd)
          cneuf(n+1)=(rn/cx)-(1.0e0_knd/cneuf(n))
60        continue
        cneuf(limn+1)=0.0e0_knd
c
c  now calculate the characteristic and exponent for each function
c  by forward multiplication using the Neumann function ratios:
        ineue(1)=int(log10(abs(cy0)))
        cneun(1)=cy0*10.0e0_knd**(-ineue(1))
        ineue(2)=int(log10(abs(cy1)))
        cneun(2)=cy1*10.0e0_knd**(-ineue(2))
          do 70 n=3,numlp
          cneun(n)=cneun(n-1)*cneuf(n-1)
          ineue(n)=int(log10(abs(cneun(n))))
          cneun(n)=cneun(n)*10.0e0_knd**(-ineue(n))
          ineue(n)=ineue(n)+ineue(n-1)
70        continue
c
c  calculate the ratios of the first derivatives of successive
c  Neumann functions using corresponding function ratios
c  calculate the ratios of the first derivatives to the
c  corresponding cylindrical Neumann functions
          do 80 n=1,limn
          rn=real(n-1,knd)
          cneudf(n)=(cx-(rn+1.0e0_knd)*cneuf(n))/(rn-cx*cneuf(n))
          cneudr(n)=(rn/cx)-cneuf(n)
80        continue
c
        if(iopr.eq.1) go to 100
c
c  convert ratios of successive functions and derivatives to ratios of
c  functions and derivatives of the same parity
          do 90 n=limn,2,-1
          cneuf(n)=cneuf(n-1)*cneuf(n)
          cneudf(n)=cneudf(n-1)*cneudf(n)
90        continue
100     continue
        return
        end
c
c
        subroutine cneuk (cm,x,limn,maxn,numlp,maxlp,ndec,nex,limbes,pi,
     1                    gamma,ngau,wr,xr,iopr,cneuf,cneun,ineue,
     2                    cneudf,cneudr)
c
c  purpose    : to calculate same parity ratios of modified cylindrical
c               Bessel functions of the third kind and the same
c               parity ratios of their first derivatives for
c               given cm and x.
c               to calculate the characteristics and exponents of the
c               modified Bessel functions of the third kind.
c               to calculate the ratios of the first derivatives to
c               corresponding modified cylindrical Bessel functions
c               of the third kind.
c
c  parameters :
c
c     input     cm     : magnitude of c, c being pure imaginary here
c               x      : x
c               limn   : the number of modified cylindrical Bessel
c                        function ratios calculated for given lnum,
c                        c, and ndec
c               maxn   : size of cneuf, cneudf and cneudr vectors
c               numlp  : the number of values of scale factors
c                        that are calculated
c               maxlp  : dimension of vectors of scale factors
c               ndec   : number of decimal digits for real(knd)
c               nex    : maximum exponent for real(knd)
c               limbes : 2*c*x+2*ndec
c               pi     : pi
c               gamma  : gamma
c               ngau   : order of the Gaussian quadrature
c               wr     : vector of ngau weighting factors for the
c                        Gaussian quadrature 
c               xr     : vector of ngau coordinates for the Gaussian
c                        quadrature
c               iopr   : if iopr = 2 convert ratios of successive 
c                        functions to ratios of successive functions of
c                        the same parity; also for the derivatives.
c
c     output  : cneuf  : same parity ratios of cylindrical Bessel
c                        functions
c               cneun  : characteristic for the cylindrical
c                        Bessel functions
c               ineue  : exponent for the cylindrical
c                        Bessel functions
c               cneudf : same parity ratios of first derivatives
c                        of cylindrical Bessel functions
c               cneudr : ratios of first derivatives to the
c                        corresponding cylindrical Bessel functions
c
        use param
c
c  real(knd) scalars and vectors
        real(knd) aj,besi0,bj,cm,cbexp,ck0,ck1,cx,dec,gamma,pi,rn,sumi,
     1            term,term1,term2,x,z
        real(knd) cbesf(limbes),cneuf(maxn),cneun(maxlp),cneudf(maxn),
     1            cneudr(maxn),wr(ngau),xr(ngau)
c
c  integer vectors
        dimension ineue(maxlp)
c
        cx=cm*x
        dec=10.0e0_knd**(-ndec-2)
c  calculate ratios of successive modified cylindrical Bessel functions
c  of the first kind
        cbesf(limbes)=cx/(2.0e0_knd*real(limbes,knd))
          do 10 n=limbes-1,1,-1
          rn=real(n+n,knd)
          cbesf(n)=1.0e0_knd/(rn/cx+cbesf(n+1))
10        continue
c
c  calculate the zeroth order function I(n=0,cm*x) divided
c  by exp(cx) by use of the ratios calculated above and the
c  formula
c
c       exp(cm*x) = I(0,cm*x) + 2I(2,cm*x) + 2I(4,cm*x) + ..., or
c
c       I(0,cm*x) = exp(cm*x)/[1 + 2(I(2,cm*x)/I(0,cm*x)) +
c                  2(I(4,cm*x)/I(0,cm*x)) + ...]
        term=2.0e0_knd
        sumi=1.0e0_knd
          do 20 j=1,limbes
          term=term*cbesf(j)
          sumi=sumi+term
          if(abs(term/sumi).lt.dec) go to 30
20        continue
30      continue
        besi0=1.0e0_knd/sumi
c
c  when cx <= 1 calculate ck0 = K(n=0,cm*x) times exp(cmx) by
c  using I(n=0,cm*x) and I ratios
        if(cx.gt.1.0e0_knd) go to 60
        term=besi0*cbesf(1)*cbesf(2)
        sumi=term
          do 40 j=2,limbes/2
          term=term*cbesf(j+j-1)*cbesf(j+j)
          sumi=sumi+term/real(j,knd)
          if(abs(term/sumi).lt.dec) go to 50
40        continue
50      ck0=(-(log(cx/2.0e0_knd)+gamma)*besi0+2.0e0_knd*sumi)
        nsub=-int(log10(abs(ck0/(ck0-2.0e0_knd*sumi))))
        if(nsub.lt.0) nsub=0
        if(nsub.gt.ndec) nsub=ndec
        nacck0=ndec-nsub-2
        ck0=ck0*exp(2.0e0_knd*cx)
        go to 100
60      if(cx.ge.40.0e0_knd) go to 70
c
c  when 40 > cm*x > 1 calculate ck0 = K(n=0,cm*x) times exp(cm*x)
c  by numerical integration (via gaussian quadrature) of its
c  integral representation
        call k0int (cx,ndec,ngau,wr,xr,ck0)
        go to 100
c
c  when cm*x >= 40 calculate ck0 = K(n=0,cm*x) times exp(cm*x) by
c  use of an asymptotic series
70      sumi=1.0e0_knd
        term1=1.0e0_knd
        z=8.0e0_knd*cx
          do 80 j=1,100
          aj=real(j,knd)
          bj=aj+aj-1.0e0_knd
          term2=-term1*bj*bj/(aj*z)
          if(abs(term2).ge.abs(term1)) go to 90
          sumi=sumi+term2
          if(abs(term2/sumi).lt.dec) go to 90
          term1=term2
80        continue
90      continue
        ck0=sqrt(pi/(2.0e0_knd*cx))*sumi
c
c  calculate ck1 = K(n=1,cm*x) times exp(cm*x) by use of the Wronskian
100     ck1=-ck0*cbesf(1)+1.0e0_knd/(cx*besi0)
c
c  now compute ratios of successive modified cylindrical Neumann
c  functions using forward recursion
c
c       K(n+1,cm*x)/K(n,cm*x)=-(2*n)/(cm*x)+1/(K(n,cm*x)/K(n-1,cm*x))
        cneuf(1)=ck1/ck0
        limit=max(numlp,limn)
          do 110 n=1,limit-1
          rn=real(n+n,knd)
          cneuf(n+1)=(rn/cx)+(1.0e0_knd/cneuf(n))
110       continue
        cneuf(limit+1)=0.0e0_knd
c
c  now calculate the characteristic and exponent for each function
c  by forward multiplication using the Neumann function ratios:
        cbexp=cx*log10(exp(1.0e0_knd))
        if(cbexp.lt.real((nex-10)*(nex-10),knd)) go to 120
        ncbexp=int(cbexp)
        cneun(1)=ck0*10.0e0_knd**(-cbexp+ncbexp)
        iterm=int(log10(abs(cneun(1))))
        cneun(1)=cneun(1)*10.0e0_knd**(-iterm)
        ineue(1)=iterm-ncbexp
        go to 150
120     numsca=cbexp/(nex-10)+1
        isca=int(cbexp/real(numsca,knd))
        cneun(1)=exp(-cx/(real(numsca,knd)))
        cneun(1)=cneun(1)*10.0e0_knd**(isca)
130     cneun(1)=ck0*(cneun(1)**(numsca))
        iterm=int(log10(abs(cneun(1))))
        cneun(1)=cneun(1)*10.0e0_knd**(-iterm)
140     ineue(1)=iterm-isca*numsca
150     cneun(2)=cneun(1)*cneuf(1)
        ineue(2)=int(log10(abs(cneun(2))))
        cneun(2)=cneun(2)*10.0e0_knd**(-ineue(2))
        ineue(2)=ineue(2)+ineue(1)
          do 160 n=3,numlp
          cneun(n)=cneun(n-1)*cneuf(n-1)
          ineue(n)=int(log10(abs(cneun(n))))
          cneun(n)=cneun(n)*10.0e0_knd**(-ineue(n))
          ineue(n)=ineue(n)+ineue(n-1)
160       continue
c
c  calculate the ratios of the first derivatives of successive
c  modified Bessel functions using corresponding function ratios
c  calculate the ratios of the first derivative to the corresponding
c  cylindrical Neumann function
          do 170 n=1,limn
          rn=real(n-1,knd)
          cneudf(n)=-(cx+(rn+1.0e0_knd)*cneuf(n))/(rn-cx*cneuf(n))
          cneudr(n)=(rn/cx)-cneuf(n)
170       continue
c
        if(iopr.eq.1) go to 190
c
c  convert ratios of successive functions and derivatives to ratios of
c  functions and derivatives of the same parity
          do 180 n=limn,2,-1
          cneuf(n)=cneuf(n-1)*cneuf(n)
          cneudf(n)=cneudf(n-1)*cneudf(n)
180       continue
190     continue
        return
        end
c
c
	subroutine gauss (n,ndec,x,w)
c
c  purpose:     To evaluate the coordinates and weighting factors
c               for an nth order Gaussian quadrature
c
c  parameters:
c
c     input:    n   : order of quadrature
c               ndec: precision for real(knd)
c
c     output:   x   : coordinate values for quadrature
c               w   : weighting factors
c
        use param
c
c  real(knd) scalars and arrays
        real(knd) delta,der,pi,s,t,test,u,v,z
        real(knd) x(n),w(n)
c
        test=10.0e0_knd**(-ndec)
        imax=(n+1)/2
        pi=acos(-1.0e0_knd)
          do 40 i=1,imax
	  z=cos(pi*(i-0.25e0_knd)/(n+0.5e0_knd))
            do 20 j=1,30
            u=0.0e0_knd
	    v=1.0e0_knd
	      do 10 k=1,n
	      t=u
              u=v
	      v=((k+k-1)*z*u-(k-1)*t)/k
10   	      continue
            s=z*z-1.0e0_knd
	    der=n*(z*v-u)/s
	    delta=-v/der-0.5e0_knd*v*v*((n*n*s-n*z*z-n)*v+
     1             2.0e0_knd*n*z*u)/(der*der*der*s*s)
            z=z+delta
	    if(abs(delta/z).lt.test) go to 30
20          continue
30        continue
	  x(i)=-z
	  x(n+1-i)=z
	  w(i)=2.0e0_knd/((1.0e0_knd-z*z)*der*der)
	  w(n+1-i)=w(i)
40	  continue
	return
	end
c
c
        subroutine k0int (cx,ndec,ngau,wr,xr,ck0)
c
c  purpose:     To evaluate the product of the modified Bessel function
c               K of order zero with argument cx times the factor
c               exp(cx)
c               Gaussian quadrature is used to evaluate the integral
c               expression for K given in 9.6.24 in the Handbook of
c               Mathematical Functions
c
c  parameters:
c
c     input:    cx   : argument of K of zero order
c               ndec : precision for real(knd)
c               ngau : order of gaussain quadrature
c               wr   : weighting factors for Gaussian quadrature
c               xr   : coordinates for Gaussian quadrature
c
c     output:   ck0  : exp(cx)*K(cx)_sub 0
c 
        use param
c
c  real(knd) scalars and vectors
        real(knd) arg,ck0,coef,cx,t,tupp
        real(knd) wr(ngau),xr(ngau)
c
        tupp=log(2.0e0_knd+4.61e0_knd*(ndec+2)/cx)
        coef=tupp/2.0e0_knd
        ck0=0.0e0_knd
          do 10 i=1,ngau
          t=coef*(1.0e0_knd+xr(i))
          arg=exp(cx*(1.0e0_knd-cosh(t)))
          ck0=ck0+wr(i)*arg
10        continue
        ck0=coef*ck0
        return
        end

