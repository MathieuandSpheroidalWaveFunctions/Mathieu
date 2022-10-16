                          Mathieu

  Matfcn is available as both a subroutine version provided as
  the module Mathieu and a stand alone version. It was originally
  developed by arnie lee van buren and jeffrey boisvert and has
  been improved several times since then.

     Table of Contents
  1. Purpose
  2. Introduction
  3. Input and Output
  4. Accuracy of results using real*8 arithmetic
  5. Accuracy of results using real*16 arithmetic
  6. Obtaining the A and B expansion coefficients
  7. Obtaining the eigenvalues


  1. Purpose

  To calculate the Mathieu radial functions and their first
  derivatives for a range of lnum orders from l = 0 to l = lnum - 1
  and for a given value of q [or c=sqrt(2*q)] and the radial
  coordinate r. The radial coordinate can either be the traditional
  coordinate z or the spheroidal-like radial coordinate xi = cosh(z).
  To calculate the Mathieu angular functions and their first derivatives
  for orders l = 0 to l = lnum - 1 for a given value of q (or c) and
  for one or more angular coordinate values.

  2. Introduction

  We use the term radial function for what is normally referred to as
  a modified Mathieu function. We do this to be consistent with our
  usage of the term radial function for spheroidal wave functions.
  The use of what we call the size parameter c instead of the
  traditional q is consistent with usage for spheroidal functions
  (see, e.g., Handbook of Mathematical Functions). So is the use of
  xi = cosh(z) as the radial coordinate. We note that the
  algebraic form of Mathieu's equation can be obtained from the
  differential equation for radial prolate spheroidal functions by
  setting m = 1/2.

  The first derivative values provided for the radial functions
  are with respect to the traditional radial coordinate z that ranges
  from 0 to infinity. First derivatives with respect to the radial
  coordinate xi can be obtained by multiplying those given by matfcn
  by the factor 1/sqrt(xi*xi-1). See the comments below about using
  x1 = xi - 1 as input to matfcn to avoid subtraction errors in forming
  the expression xi*xi - 1 when xi is near unity. Here xi*xi - 1 is
  then accurately computed using x1*(x1 + 2).

  The first derivative values provided for the angular functions
  are with respect to the traditional angle coordinate.

  The methods used in matfcn are documented in the journal article:
  A. L. Van Buren and J. E. Boisvert, "Accurate calculation of the
  modified mathieu functions of integer order," Quart. Appl.
  Math. 65, 1-23 (2007). A pdf version of this article is provided
  in this repository. [Note that the argument for both cos and sin
  in the two equations in (10) is in error and should be r times
  theta.]

  Matfcn is written in free format fortran. It is designed around the
  maximum number of decimal digits ndec and the maximum exponent nex
  available in real arithmetic. Procedures used in matfcn allow for
  exponents much larger than nex since the resulting floating point
  radial function values are given as a characteristic and an integer
  exponent. Testing of matfcn shows that in real*8 arithmetic,it provides
  good results for q positive for virtually all values of xi and for c
  up to at least 10000. It also provides good results for q negative for
  virtually all values of xi down to 1.000005 and for c up to at
  least 10000i. See the discussion of accuracy below.

  Matfcn can be run in either double precision or quadruple precision
  arithmetic. The choice is set in the module param provided in the github
  repository. If this is not available, then create param as follows:

    module param
    integer, parameter :: knd = selected_real_kind(8)
    logical, parameter :: debug = .true.
    logical, parameter :: warn = .true.
    logical, parameter :: output = .true.
    end module param

  Set the value of knd in the parenthesis to either 8 for double
  precision or 16 for quadruple precision arithmetic. Using quadruple
  precision will provide higher accuracy but will considerably increase
  the run time. Some compilers require that param be compiled prior to
  rather than after mathieu. The logicals in param are described below
  in the discussion of the output files.

  Some computers may have more than 8 bytes for double precision
  data and more than 16 bytes for quadruple precision data. In this
  case just use the appropriate integers for the kind parameters in
  module param. Also change the values of kindd and kindq set in
  statement 5 below below the comments section to the number of
  bytes for double precision data and quadruple precision data,
  respectively. Larger values of kindd and kindq will lead to
  more accurate results, especially for large values of c and small
  values of x. Note that the mathematical constant gamma is given below
  to 39 decimal digits. Of course, this will be truncated to the
  correct precision for the kind parameter used. If matfcn is run
  with higher precision than this, additional digits will be necessary
  for gamma in order to ensure the higher accuracy associated with
  this higher precision.

  Some computer also use values for kind that do not correspond to the
  number of bytes used in real data, In this case just use the values
  for kind that correspond to double and quadruple precison arithmetic.
  This includes setting kindd and kindq to the proper values for double
  and quadruple precision arithmetic, respectively.

  3. Input and Output

  Following is a description of the input and output parameters in the
  call statement for the subroutine version. After that will be a
  description of the the input and output files associated with the
  stand alone version. Note that these output files, if desired, can be
  obtained when running the subroutine version. See comments about
  this below.

  A sample input and resulting output from matfcn for c real (q positive)
  is provided by the files matfcndat (text version of the input file
  matfcn.dat for the stand alone version), matfort20 (text version of the
  output file fort.20 giving the resulting radial functions) and matfort30
  (text version of the output file fort.30 giving the resulting angular
  functions).

  A sample input and resulting output from matfcn for c imagimary (q negative)
  is provided by the files matfcnidat (text version of the input file
  matfcn.dat for the stand alone version), matifort20 (text version of the
  output file fort.20 giving the resulting radial functions) and matifort30
  (text version of the output file fort.30 giving the resulting angular
  functions).

  Subroutine Version of Matfcn

    subroutine matfcn(lnum, ioprad, izxi, icq, isq, qc, r, iopang, narg, arg, &
                      mc1c, mc1e, mc1dc, mc1de, mc23c, mc23e, mc23dc, mc23de, naccrc, &
                      ms1c, ms1e, ms1dc, ms1de, ms23c, ms23e, ms23dc, ms23de, naccrs, &
                      ce, ced, se, sed, nacca)

         integer, intent (in)    ::  lnum, ioprad, izxi, icq, isq, iopang, narg
         real(knd), intent (in)  ::  qc, r, arg(narg)
         integer, intent (out)   ::  mc1e(lnum), mc1de(lnum), mc23e(lnum), mc23de(lnum), &
                                     ms1e(lnum), ms1de(lnum), ms23e(lnum), ms23de(lnum), &
                                     naccrc(lnum), naccrs(lnum), nacca(lnum, narg)
         real(knd), intent (out) ::  mc1c(lnum), mc1dc(lnum), mc23c(lnum), mc23dc(lnum), &
                                     ms1c(lnum), ms1dc(lnum), ms23c(lnum), ms23dc(lnum), &
                                     ce(lnum, narg), ced(lnum, narg), &
                                     se(lnum, narg), sed(lnum, narg)

    Input and output parameters appearing in the subroutine call
    statement are defined below:

          lnum   : number of integer values of l, given by
                   0, 1, 2, ..., lnum-1, that radial and/or
                   angular functions are desired

          ioprad : (integer)
                 : =0 if radial functions are not desired
                 : =1 if radial functions of both kinds and
                      their first derivatives are computed.
                      When q is positive, this results in
                      radial functions of the first kind Mc1 and
                      Ms1 and of the second kind Mc2 and Ms2 and
                      their first derivatives. When q is negative,
                      the radial functions of the second kind are
                      replaced by those of the third kind Mc3
                      and Ms3 and their first derivatives.

          izxi   : (integer)
                 : = 1 if the input radial coordinate is the
                       traditional radial coordinate z
                 : = 2 if the input radial coordinate is the
                       spheroidal-like radial coordinate = cosh(z)
                       (actually x1 = xi - 1 is imput to avoid
                       subtraction errors in the calculation of
                       xi*xi - 1 in matfcn when xi is near unity.
                       We avoid the subtraction error by using
                       instead x1*(x1 + 2). Note that When z is
                       input, the subtraction error is avoided in
                       computing x1 = xi - 1 by using the identity
                       x1=2*sinh(z/2)*sinh(z/2)

          icq    : integer equal to 1 when q is input and equal to 2
                   when cm, the magnitude of c, is input
                   Note that c is pure imaginary when q is negative;
                   we choose the positive square root of q in the
                   definition of c.

          isq    : integer equal to +1 when q is positive and
                   equal to -1 when q is negative

          qc     : size parameter [real(knd)]
                 : = q when icq = 1
                 : = cm, magnitude of c, when icq = 2
                   This is equal to c when isq = +1 so that q
                   is positive. When isq = -1 so that q
                   is negative, then c is equal to i times the
                   input value

          r      : radial coordinate [real(knd)]
                 : = z [real(knd)] when izxi = 1; this is the
                   traditional radial coordinate ranging from
                   zero to infinity.
                 : = xi - 1 [real(knd)] when izxi = 2; where xi is
                   the spheroidal-like radial coordinate equal to
                   cosh(z). A dummy value, e.g. 1.0d0, should be
                   entered for r if ioprad = 0

          iopang : (integer)
                 : = 0 if angular functions are not desired
                 : = 1 if angular functions of the first kind
                       are computed
                 : = 2 if angular functions of the first kind and
                       their first derivatives are computed

          narg   : number of values of the angular coordinate phi for
                   which angular functions of the first kind are
                   calculated (integer). If no angular functions
                   are desired (iopang = 0), set narg equal to unity.

          arg:     vector arg(narg) containing the values of phi in
                   degrees for which angular functions are desired
                   [real(knd)]

          mc1c   : real(knd) vectors of length lnum containing the
          mc1dc    characteristics for the cosine radial functions
                   of the first kind Mc1 and their first derivatives
                   with respect to z.
                   When q is positive, the functions are real.
                   When q is negative and the order is odd, the
                   functions are imaginary. The real values given
                   in the vectors mc1c and mc1dc must be multiplied
                   by i to obtain the characteristics for odd orders.

          mc1e   : integer vectors of length lnum containing the
          mc1de    exponents corresponding to mc1c and mc1dc.

          mc23c  : When q is positive, these are real(knd) vectors of
          mc23dc   length lnum containing the characteristics for the
                   cosine radial functions of the second kind Mc2 and
                   their first derivatives with respect to z.
                   When q is negative, these vectors instead contain
                   the characteristics for the cosine radial functions
                   of the third kind Mc3 and their first derivatives
                   with respect to z. Mc3 and its first derivative is
                   real when l is odd and imaginary when l is even.
                   In this case, the real values given in the vectors
                   mc23c and mc23dc must be multiplied by i to obtain
                   the characteristics for even orders.

          mc23e :  integer vectors of length lnum containing the
          mc23de   exponents corresponding to mc23 and mc23d

          naccrc : vector of lnum values for the estimated accuracy of the
                   cosine radial functions; obtained using the Wronskian

          ms1c   : real(knd) vectors of length lnum containing the
          ms1dc    characteristics for the sine radial functions
                   of the first kind Ms1 and their first derivatives
                   with respect to z.
                   When q is positive, the functions are real.
                   When q is negative and the order is odd, the
                   functions are imaginary. The real values given in
                   the vectors ms1c and ms1dc must be multiplied by i
                   to obtain the characteristics for odd orders.

          ms1e   : integer vectors of length lnum containing the
          ms1de    exponents corresponding to ms1c and ms1de

          ms23c  : When q is positive, these are real*8 vectors of
          ms23dc   length lnum containing the characteristics for the
                   sine radial functions of the second kind Ms2 and
                   their first derivatives with respect to z.
                   When q is negative, these vectors instead contain
                   the characteristics for the sine radial functions
                   of the third kind Ms3 and their first derivatives
                   with respect to z. Ms3 and its first derivative is
                   real when l is odd and imaginary when l is even.
                   In this case, the real values given in the vectors
                   ms23c and ms23dc must be multiplied by i to obtain
                   the characteristics for even orders.

          mc23e  : integer vectors of length lnum containing the
          mc23de   exponents corresponding to ms23 and ms23d

          naccrs : vector of lnum values for the estimated accuracy of the
                   sine radial functions; obtained using the Wronskian

          ce,ced : arrays ce(lnum,narg) and ced(lnum,narg) that
                   contain narg calculated angular cosine functions
                   and their first derivatives for each of the lnum
                   values of l [real(knd)]
                   For example, ce(10,1) is the angular function for
                   l = 9 and the first value of the angle phi given
                   by arg(1)

          se,sed : arrays se(lnum,narg) and sed(lnum,narg) that
                   contain narg calculated angular sine functions
                   and their first derivatives for each of the lnum
                   values of l [real(knd))
                   For example, se(10,1) is the angular function for
                   l = 9 and the first value of the angle phi given
                   by arg(1)

          nacca  : array of lnum values for the estimated accuracy of the angular
                   functions for each of the narg angles; it is equal to the
                   minimum of the estimates for the sine and cosine angular
                   functions and their first derivatives; it is based on
                   subtraction errors in their calculation

  Stand Alone Version of Matfcn

        Input Data

  Input parameters are read from unit 1 in the file matfcn.dat
  The user should create an input file with the lines of data as
  described.

       line 1:

          lnum   : number of integer values of l given by
                   0, 1, 2, ..., lnum-1 (integer)

          ioprad : (integer)
                 : = 0 if radial functions are not desired
                 : = 1 if radial functions of both kinds and
                       their first derivatives are computed.
                       When q is positive, this results in
                       radial functions of the first and second
                       kind m1 and m2. When q is negative, it
                       results in radial functions of the first
                       and third kind m1 and m3.

          iopang : (integer)
                 : = 0 if angular functions are not desired
                 : = 1 if angular functions of the first kind
                       are computed
                 : = 2 if angular functions of the first kind and
                       their first derivatives are computed

          izxi   : (integer)
                 : = 1 if the input radial coordinate is the
                       traditional radial coordinate z
                 : = 2 if the input radial coordinate is the
                       spheroidal-like radial coordinate xi = cosh(z)
                       (actually x1 = xi - 1 is imput to avoid
                       subtraction errors in the calculation of
                       xi*xi - 1 in matfcn when xi is near unity.
                       We avoid the subtraction error by using
                       instead x1*(x1 + 2). Note that When z is
                       input, the subtraction error is avoided in
                       computing x1 = xi - 1 by using the identity
                       x1=2*sinh(z/2)*sinh(z/2)

          icq    : integer equal to 1 when q is input and equal to 2
                   when cm, the magnitude of c, is input
                   Note that c is pure imaginary when q is negative;
                   we choose the positive square root of q in the
                   definition of c.

          isq    : integer equal to +1 when q is positive and
                   equal to -1 when q is negative. This is used
                   to determine q when cm is input.

       line 2:

          q or cm: input size parameter [real(knd)]
                 : q when icq = 1
                 : cm, magnitude of c, when icq=2
                   This is equal to c when isq = +1 so that q
                   is positive. When isq = -1 so that q
                   is negative, then c is equal to i times the
                   input value.

       line 3:

          r      : input radial coordinate [real(knd)]
                   for izxi = 1: equal to the traditional radial
                   coordinate z ranging from zero to infinity.
                   for izxi = 2: equal to xi - 1 where xi is the
                   spheroidal-like radial coordinate equal to cosh(z).
                   A dummy value, e.g. 1.0d0, should be entered for r
                   if no radial functions are desired (ioprad=0)

       line 4:

          arg1   : first value for the angle coordinate (in degrees)
                   [real(knd)]; only read if iopang not equal to 0

          darg   : increment used to calculate additional desired
                   arguments for angular functions. [real(knd)]

          narg   : number of desired angle arguments. (integer)
                   enter the value 1 for narg if iopang = 0


         Output files

   These output files are also available using the subroutine version
   of matfcn. Generation of each of the files is controlled by a logical
   specified in the module param. False suppresses the output file and true 
   enables it. The logical debug controls fort.30 and fort.40, the logical
   output controls fort.20 and fort.30 and warn controls fort.60.

   fort.20

     This file contains values for all radial functions that have
     been calculated. The first line contains the values for either
     z or xi, depending on whether z or xi - 1 is input and either c
     (the factor i = sqrt(-1) is suppressed from c when q is negative)
     or q, depending on whether c or q is input, formatted as follows
     (see statements 20 and 30 in the first subroutine main):

                z or xi : e24.15 for real*8; e40.31 for real*16
                c or q  : e24.15 for real*8; e40.31 for real*16

     Each subsequent line in matfcn.rad contains radial functions
     for given values of l, beginning with 0 and continuing to l =
     lnum - 1. First is given the value of l, followed on the
     same line first by the cosine radial functions and then by a
     two digit estimate of the number of accurate digits in the
     preceeding radial functions. For l unequal to zero, the line
     of cosine radial functions is followed by a line of the
     corresponding sine radial functions for that value of l
     together with a two digit estimate of their accuracy.
     the estimate of accuracy is determined by comparing the
     theoretical Wronskian to the value computed using the radial
     functions and their first derivatives. The output and
     corresponding format for each line is as follows
     (see statements 520, 530, 550, and 560 in subroutine main)

        for ioprad = 1

               l      : value for l (i5)
               mc1c/  : characteristic of the cosine/sine radial
               ms1c     function of first kind (f17.14). when q is
                        negative (c is imaginary) and l is odd, mc1c
                        and ms1c are pure imaginary and the factor i
                        is suppressed
               imc1e/ : exponent of the cosine/sine radial function of
               ims1e    first kind (i6 for q positive; i7 for q
                        negative)
               mc1dc/ : characteristic of the first derivative of the
               ms1dc    radial function of first kind (f17.14). when q
                        is negative and l is odd, mc1dc and ms1dc are
                        pure imaginary and the factor i is suppressed
               imc1e/ : exponent of the first derivative of the
               ims1e    radial function of first kind (i6 for q
                        positive; i8 for q negative)
               mc23c/ : characteristic of the radial function
               ms23c    of second (or third) kind (f17.14). when q
                        is negative and l is even, mc3c and ms3c are
                        pure imaginary and the factor i is suppressed
               imc23e/: exponent of the radial function of the
               ims23e   second (or third) kind (i6 for q positive;
                        i8 for q negative)
               mc23dc/: characteristic of the first derivative of the
               ms23dc   radial function of the second (or third) kind
                        (f17.14). when q is negative and l is even,
                        mc3dc and ms3dc are pure imaginary and the
                        factor i is suppressed
               imc23de/:exponent of the first derivative of the
               ims23de  radial function of the second (or third)
                        kind (i6 for q positive; i8 for q negative)
               naccrc/ :accuracy: equal to the number of decimal digits
               naccrs   of agreement between the theoretical Wronskian
                        and the calculated Wronskian (i2)

     When one of the two products of radial functions used in forming
     the Wronskian is significantly smaller in magnitude than the other
     product, the Wronskian accuracy naccr might overestimate the
     accuracy of one or both of the two functions in the smaller
     product. They can be less accurate than naccr by as many digits
     as obtained by integer truncation of the logarithm to the base 10
     of the ratio of the magnitude of the larger product to the smaller
     product. The file fort.40 (described below) will show the
     subtraction errors involved in the calculation of the functions
     and identify inaccuracies in their values.

   fort.30

     This file contains values for all angular functions that have
     been calculated. The first line contains the value for either c
     (i is suppressed when q is negative) or q, depending on which
     one is input, formatted as follows (see statements 40 and 45
     in subroutine main):

                c or q : e24.15 in real*8; e40.31 for real*16

     The second line contains the value for the first l (=0),
     formatted as follows (see statement 90 in subroutine main):

                l      : i5

     This is followed by a series of narg lines. each line contains
     a desired value of angle followed by the corresponding angular
     functions and accuracy. Specific output and format for each line
     is as follows:

        for iopang = 1:

               arg    : for ioparg = 0, angle in degrees (f20.14; see
                        statements 620 and 630 in subroutine main)
               ce     : value of the angular function ce (e24.15;
                        see statements 620 and 630 in subroutine main)
               se     : value of the angular function se (e24.15;
                        see statements 620 and 630 in subroutine main).
                        for l = 0 the value for se is set equal to zero.

        for iopang = 2, each line also includes:
               ced    : value of the first derivative of the
                        angular function ce, listed immediately
                        following the value for ce: (e24.15;
                        see statement 630 in subroutine main)
               sed    : value of the first derivative of the
                        angular function se, listed immediately
                        following the value for se: (e24.15);
                        see statement 630 in subroutine main).
                        for l = 0 the value for sed is set equal to
                        zero.

        for iopang = 1 or 2:
      nacca           : accuracy: smaller of naccc and naccs, where
                        naccc is the number of decimal digits of
                        accuracy estimated for the angular function
                        ce and its first derivative when iopang = 2
                        and where naccs is the corresponding estimate
                        for se. naccc and naccs are conservative
                        estimates given by ndec -2 - the larger of
                        the subtraction error encountered in the
                        series calculation of the corresponding
                        angular function and its first derivative.
                        When nacca is equal to 0, the corresponding
                        angular functions and their first derivatives
                        are set equal to zero. naccc and naccs are
                        written with format i2; see statements 620
                        and 630 in subroutine main)

   fort.40 and fort.50

     These files are diagnostic files. They contain information about
     the specific techniques used and the numbers of terms required
     for the radial function and angular function calculations,
     respecively. They are annotated and should be self-explanatory.

   fort.60

     This file may be of interest to the user. It is recommended that
     the user utilize this file, especially when using matfcn for very
     large imaginary values of c with x near unity. Whenever the
     estimated accuracy for the radial functions falls below a
     designated integer value, the associated values of x, c (and the
     corresponding value for q), and l are written to fort.60. This
     integer is currently set equal to 6 in the write statement for
     fort.60 found just before the line numbered 600 below in
     subroutine main. The reason for choosing 6 is that it is expected
     that 6 accurate decimal digits are sufficient for most
     applications. The integer can be changed to any desired value in
     the write(60,*).

     In the unlikely case that the eigenvalue routine fails to properly
     converge, the value of l for which this occurs will be written
     to fort.60. Note that this has never been observed with matfcn.


  4. Accuracy of Results Using Real*8 Arithmetic

  When c is real (q positive) and xi is not zero, matfcn provides
  values for the radial functions of both the first and second kinds
  that are usually accurate to at least 10 decimal digits. Matfcn
  estimates the accuracy as the nummber of decimal digits of agreement
  between the theoretical Wronskian and its value computed using the
  resulting radial function values. We tested matfcn for a selection
  of xi values ranging from 1.00000000001 to 100 and for a range of c
  values up to 1000. We also tested all of the xi values up to 10 for
  c values up to 10000. Lnum was chosen to be the larger of 1000 or an
  order where the radial functions of the first kind had magnitudes
  well below 10 to the -300 power. The resulting accuracies were 10 or
  more digits with c up to 100. Accuracies were at least 8 digits for c
  = 200 with only three occurrences of 8 digits. For c = 500, there
  were three 7 digit results. For larger c there were a few 7 digit
  and two 6 digit results. Matfcn is expected to provide useful results
  for xi outside the test range and for lnum much higher than tested.
  It is always possible to choose xi so close to a root that lower
  accuracy is obtained. The corresponding radial function would then
  have a somewhat smaller magnitude than neighboring function values
  and likely would not adversely affect solutions to problems involving
  these functions.

  When c is real and xi is equal to unity, the sine radial functions
  of the first kind and the first derivative of the cosine radial
  function of the first kind are equal to zero. The other radial
  functions of the first kind are very accurate unless extremely close
  to a root. When c is large and the order l is somewhat less than
  2c/pi, both the cosine radial functions of the second kind and the
  first derivatives of the sine radial functions of the second kind
  tend to have a very small magnitude. For a given large value of c,
  the magnitude decreases as the order decreases toward zero. For a
  given order, the magnitude decreases as c increases. These functions
  are inaccurate to the degree that they are small. Thus they should
  not reduce the accuracy of solutions to problems where the radial
  functions of the second kind appear in a pair with the functions of
  the first kind. Here the accurate radial functions of the first kind
  are larger and dominate. Note that the Wronskian accuracy check
  overestimates the accuracy here.

  When c is imaginary (q negative) and xi is not equal to zero, matfcn
  provides values for the radial functions of the first and third kinds
  that are usually accurate to 10 or more digits. We tested matfcn
  using a set of xi values ranging from 1.000000000001 to 100 and a set
  of c values up to 10000i. We tested the full set of xi values for c
  up to 10i. We tested the xi values ranging from 1.000005 to 100 for c
  up to 5000i and the xi values ranging from 1.000005 to 10 for c up to
  10000i. Lnum was chosen to be the larger of 1000 or a value large
  enough that the radial functions of the first kind have magnitudes
  well below 10 to the -300 power. The resulting accuracies were 10 or
  more digits for c up to 50i and 9 or more digits for c up to 500i.
  For higher values of c, accuracies were 8 or more digits except for a
  few 7 digit results and one 6 digit result.

  For c larger than about 10i, matfcn uses the Bessel function (K)
  expansion to compute radial functions of the third kind for orders
  less than 2*cm/pi. The Bessel product expansion suffers too much
  subtraction error here to provide accurate results. The number of
  terms maxk required for convergence of the K expansion to the
  desired accuracy and the time to perform the calculations increases
  somewhat linearly with 1/(xi-1). For xi = 1.000001 and c = 10000i,
  storing the vectors of expansion functions and coefficients for the
  calculations can require as much as 1 GB of memory. I set the maximum
  value of maxk to 10000000 (maxterm) for my testing. The user can
  increase this in statement number 20 below if they desire smaller xi
  values for c greater than 10i and have sufficient computer memory.
  Note that the necessary memory doubles for real*16 arithmetic and
  that the Bessel function expansion can be used to higher values
  of c before the K expansion is required.

  When c is imaginary and xi is equal to unity, the cosine radial
  function of the first kind and the first derivative of the sine
  radial function are equal to zero. The Bessel function expansion
  used for radial functions of the third kind for orders below about
  2(cm)/pi is useless for xi = 0. Unfortunately the Bessel product
  expansion that must be used here suffers subtraction errors that
  are too large to allow accurate values for low order radial functions
  of the third kind when c is somewhat greater than 20i. Note
  that when the function accuracy is such that the Wronskian accuracy
  is less than the desired minimum accuracy, matfcn uses the
  Wronskian to compute an accurate value for the sine radial
  function of the third kind and the first derivative of the cosine
  radial function of the third kind. Of course, the cosine radial
  function and the first derivative of the sine radial function
  likely have an accuracy less than the desired minimum accuracy and
  in fact are set equal to zero when the Wronskian aaccuracy is zero.
  Radial function values for orders near and above 2(cm)/pi are
  usually highly accuracy for all values of c.

  The calculated angular functions are highly accurate except for
  orders less than about 2(cm)/pi when c is large in magnitude.
  Function values tend to be less accurate the larger the magnitude of
  c, the smaller the order and the closer the angle is to zero degrees.
  However, the loss in accuracy is due to subtraction error in the
  calculations and is accompanied by a proportional decrease in the
  magnitude of the angular function value relative to the corresponding
  trigonometric function values. This decrease in magnitude almost
  always results in a corresponding reduction in the magnitude of their
  contribution to the solution of physical problems involving Mathieu
  functions. Note that the angular functions have the same norm pi as
  the corresponding cosine and sine functions.

  5. Accuracy of Results Using Real*16 Arithmetic

  If the user desires more accuracy than is provided using real*8
  arithmetic, matfcn can be run using real*16 arithmetic. However,
  the execution time increases considerably. Matfcn computes results
  to full precision except for the radial functions of the third kind
  m3 for orders less than 2*cm/pi. These are computed to a precision
  based on the parameter minacc. Minacc is set below after these
  comments to be equal to 15 digits for real*16 arithmetic. It can be
  increased to provide greater precision and accuracy for lower order
  values of m3. Although when c is imaginary and large in magnitude and
  xi is very near unity, the accuracy may not increase much when using
  the expansion in K functions. Here, many more terms than used in
  matfcn would be required to obtain significantly higher accuracy.
  It is advised that one not change minacc from its value of 10 given
  below when using real*8 arithmetic.

  If desired, output data can be written to one or more of the
  following five files: fort.20, fort.30, fort.40, fort.50 and fort.60.
  The write and corresponding format statements are commented out
  below to allow users to avoid these files. The user can obtain any
  desired output file, say fort.20, by searching matfcn for write(20,
  and the associated format statements and removing the c in column 1
  of the lines of code found. Some compilers also require that the file
  be opened. In that case remove the c in the line of code below
  opening the file. Users running matfcn outside the range over which
  it was tested, especially when the size parameter c has a magnitude
  larger than 10000 may want to consider using fort.60. Matfcn writes
  to fort.60 the values for xi, c, m, l, and estimated accuracy when
  the estimated accuracy falls below a specified number of integer
  digits. See comments below about this file.

  6. Obtaining the expansion coefficients A and B

  The user may be interested in the expansion coefficients A and B for
  each order l. Both are calculated as vectors containing ratios of
  successive coefficients. The vector enra(j) contains the ratios
  A(sub 2j)/A(sub 2j-2) for l even and A(sub 2j+1)/A(sub 2j-1) for l
  odd. Enra is computed in convera and returned to subroutine main.
  Similarly, the vector enrb(j) contains the ratios B(sub 2j)/
  B(sub 2j-2) for even l and B(sub 2j+1)/B(sub 2j-1) for l odd. Enrb is
  computed in converb and also returned to subroutine main. Note that
  B(sub 0) is not used so the ratio enrb(1) for l even is also not
  used.

  The subroutines convera and converb return coefficients appropriate
  for q positive. If q is negative, the A and B coeffficients are the
  same as those for q positive when l is even. However, when l is odd,
  the A coefficients for q negative are equal to the B coefficients for
  q positive. And the B coefficients for q negative are equal to the A
  coefficients for q positive. Thus the coefficient ratios enra and
  enrb returned by convera and converb must be interchanged when q is
  negative.

  Both enra and enrb have the dimension maxd, which is the
  estimate of the largest number of coefficients that will be required
  to compute the radial and/or angular functions for all desired lnum
  orders. The number of coeffficient ratios computed for each l is
  given by limd/2, where limd tends to increase with increasing l. It
  can be very large when xi is very close to unity and the traditional
  expansion in modified Bessel function of the second kind K is used to
  compute the radial functions m3. If enra and enrb are to be used in
  calculations involving just the the angular functions, then only
  limdang/2 coefficients are needed for each l. The value of limdang is
  determined for each l below between statements 90 and 95 in
  subroutine main. Note that limdang is only determined when matfcn is
  run with iopang equal to 1 or 2.

  Matfcn computes the coefficient A(sub l) for order l in subroutine
  cese where the angular functions are computed and returns it to
  subroutine main as asubl. The A coefficient with the next higher
  index A(sub l+2) is obtained by multiplying A(sub l) by
  enra(sub (l+2)/2) for l even and by enra(sub (l+1)/2) for l odd.
  Higher index cofficients are then obtained in order by multiplying by
  successive ratios enra. The A coefficient A(sub l-2) is obtained by
  dividing A(sub l) by the ratio enra(sub l/2) for even l and by
  enra(sub (l-1)/2) for l odd. Lower order coefficients are then
  obtained in turn by successive division by lower index ratios enra.
  The B coefficients are obtained in a similar fashion starting with
  B(sub l) computed in subroutine cese and returned to main as bsubl.
  Note that the values for asubl and bsubl are correct as given for
  both q positive and q negative. No interchange is needed when q is
  negative and l is odd.

  7. Obtaining the eigenvalues

  The eigenvalues for the expansion coefficients A are computed in
  subroutine convera and returned to main where they are stored in the
  vector eiga(l+1). Similarly the eigenvalues for B are computed in
  subroutine converb and returned to main and stored in the vector
  eigb(l+1). Note that there is no B eigenvalue for l = 0 so eigb(1)
  is set to zero here. Note also that for q negative the eigenvalues
  for A and B switch places when l is odd; thus the A eigenvalue is
  computed in converb and the B eigenvalue is computed in convera. The
  vectors eiga and eigb for q negative given below correctly reflect
  this interchange so no further change needs to be made.
