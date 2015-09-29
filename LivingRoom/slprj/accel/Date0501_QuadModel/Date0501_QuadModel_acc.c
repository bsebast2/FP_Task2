#include "__cf_Date0501_QuadModel.h"
#include <math.h>
#include "Date0501_QuadModel_acc.h"
#include "Date0501_QuadModel_acc_private.h"
#include <stdio.h>
#include "simstruc.h"
#include "fixedpoint.h"
#define CodeFormat S-Function
#define AccDefine1 Accelerator_S-Function
#ifndef __RTW_UTFREE__  
extern void * utMalloc ( size_t ) ; extern void utFree ( void * ) ;
#endif
boolean_T Date0501_QuadModel_acc_rt_TDelayUpdateTailOrGrowBuf ( int_T *
bufSzPtr , int_T * tailPtr , int_T * headPtr , int_T * lastPtr , real_T
tMinusDelay , real_T * * tBufPtr , real_T * * uBufPtr , real_T * * xBufPtr ,
boolean_T isfixedbuf , boolean_T istransportdelay , int_T * maxNewBufSzPtr )
{ int_T testIdx ; int_T tail = * tailPtr ; int_T bufSz = * bufSzPtr ; real_T
* tBuf = * tBufPtr ; real_T * xBuf = ( NULL ) ; int_T numBuffer = 2 ; if (
istransportdelay ) { numBuffer = 3 ; xBuf = * xBufPtr ; } testIdx = ( tail <
( bufSz - 1 ) ) ? ( tail + 1 ) : 0 ; if ( ( tMinusDelay <= tBuf [ testIdx ] )
&& ! isfixedbuf ) { int_T j ; real_T * tempT ; real_T * tempU ; real_T *
tempX = ( NULL ) ; real_T * uBuf = * uBufPtr ; int_T newBufSz = bufSz + 1024
; if ( newBufSz > * maxNewBufSzPtr ) { * maxNewBufSzPtr = newBufSz ; } tempU
= ( real_T * ) utMalloc ( numBuffer * newBufSz * sizeof ( real_T ) ) ; if (
tempU == ( NULL ) ) { return ( false ) ; } tempT = tempU + newBufSz ; if (
istransportdelay ) tempX = tempT + newBufSz ; for ( j = tail ; j < bufSz ; j
++ ) { tempT [ j - tail ] = tBuf [ j ] ; tempU [ j - tail ] = uBuf [ j ] ; if
( istransportdelay ) tempX [ j - tail ] = xBuf [ j ] ; } for ( j = 0 ; j <
tail ; j ++ ) { tempT [ j + bufSz - tail ] = tBuf [ j ] ; tempU [ j + bufSz -
tail ] = uBuf [ j ] ; if ( istransportdelay ) tempX [ j + bufSz - tail ] =
xBuf [ j ] ; } if ( * lastPtr > tail ) { * lastPtr -= tail ; } else { *
lastPtr += ( bufSz - tail ) ; } * tailPtr = 0 ; * headPtr = bufSz ; utFree (
uBuf ) ; * bufSzPtr = newBufSz ; * tBufPtr = tempT ; * uBufPtr = tempU ; if (
istransportdelay ) * xBufPtr = tempX ; } else { * tailPtr = testIdx ; }
return ( true ) ; } real_T Date0501_QuadModel_acc_rt_TDelayInterpolate (
real_T tMinusDelay , real_T tStart , real_T * tBuf , real_T * uBuf , int_T
bufSz , int_T * lastIdx , int_T oldestIdx , int_T newIdx , real_T initOutput
, boolean_T discrete , boolean_T minorStepAndTAtLastMajorOutput ) { int_T i ;
real_T yout , t1 , t2 , u1 , u2 ; if ( ( newIdx == 0 ) && ( oldestIdx == 0 )
&& ( tMinusDelay > tStart ) ) return initOutput ; if ( tMinusDelay <= tStart
) return initOutput ; if ( ( tMinusDelay <= tBuf [ oldestIdx ] ) ) { if (
discrete ) { return ( uBuf [ oldestIdx ] ) ; } else { int_T tempIdx =
oldestIdx + 1 ; if ( oldestIdx == bufSz - 1 ) tempIdx = 0 ; t1 = tBuf [
oldestIdx ] ; t2 = tBuf [ tempIdx ] ; u1 = uBuf [ oldestIdx ] ; u2 = uBuf [
tempIdx ] ; if ( t2 == t1 ) { if ( tMinusDelay >= t2 ) { yout = u2 ; } else {
yout = u1 ; } } else { real_T f1 = ( t2 - tMinusDelay ) / ( t2 - t1 ) ;
real_T f2 = 1.0 - f1 ; yout = f1 * u1 + f2 * u2 ; } return yout ; } } if (
minorStepAndTAtLastMajorOutput ) { if ( newIdx != 0 ) { if ( * lastIdx ==
newIdx ) { ( * lastIdx ) -- ; } newIdx -- ; } else { if ( * lastIdx == newIdx
) { * lastIdx = bufSz - 1 ; } newIdx = bufSz - 1 ; } } i = * lastIdx ; if (
tBuf [ i ] < tMinusDelay ) { while ( tBuf [ i ] < tMinusDelay ) { if ( i ==
newIdx ) break ; i = ( i < ( bufSz - 1 ) ) ? ( i + 1 ) : 0 ; } } else { while
( tBuf [ i ] >= tMinusDelay ) { i = ( i > 0 ) ? i - 1 : ( bufSz - 1 ) ; } i =
( i < ( bufSz - 1 ) ) ? ( i + 1 ) : 0 ; } * lastIdx = i ; if ( discrete ) {
double tempEps = ( DBL_EPSILON ) * 128.0 ; double localEps = tempEps *
muDoubleScalarAbs ( tBuf [ i ] ) ; if ( tempEps > localEps ) { localEps =
tempEps ; } localEps = localEps / 2.0 ; if ( tMinusDelay >= ( tBuf [ i ] -
localEps ) ) { yout = uBuf [ i ] ; } else { if ( i == 0 ) { yout = uBuf [
bufSz - 1 ] ; } else { yout = uBuf [ i - 1 ] ; } } } else { if ( i == 0 ) {
t1 = tBuf [ bufSz - 1 ] ; u1 = uBuf [ bufSz - 1 ] ; } else { t1 = tBuf [ i -
1 ] ; u1 = uBuf [ i - 1 ] ; } t2 = tBuf [ i ] ; u2 = uBuf [ i ] ; if ( t2 ==
t1 ) { if ( tMinusDelay >= t2 ) { yout = u2 ; } else { yout = u1 ; } } else {
real_T f1 = ( t2 - tMinusDelay ) / ( t2 - t1 ) ; real_T f2 = 1.0 - f1 ; yout
= f1 * u1 + f2 * u2 ; } } return ( yout ) ; } real_T
rt_urand_Upu32_Yd_f_pw_snf ( uint32_T * u ) { uint32_T lo ; uint32_T hi ; lo
= * u % 127773U * 16807U ; hi = * u / 127773U * 2836U ; if ( lo < hi ) { * u
= 2147483647U - ( hi - lo ) ; } else { * u = lo - hi ; } return ( real_T ) *
u * 4.6566128752457969E-10 ; } real_T rt_nrand_Upu32_Yd_f_pw_snf ( uint32_T *
u ) { real_T y ; real_T sr ; real_T si ; do { sr = 2.0 *
rt_urand_Upu32_Yd_f_pw_snf ( u ) - 1.0 ; si = 2.0 *
rt_urand_Upu32_Yd_f_pw_snf ( u ) - 1.0 ; si = sr * sr + si * si ; } while (
si > 1.0 ) ; y = muDoubleScalarSqrt ( - 2.0 * muDoubleScalarLog ( si ) / si )
* sr ; return y ; } static void mdlOutputs ( SimStruct * S , int_T tid ) {
real_T jq5l4n1lxp ; real_T * lastU ; real_T fegfldezgg [ 3 ] ; real_T
gg2u4if1k1 [ 3 ] ; real_T h4fomgvj5n ; real_T erv0ccpdoz ; real_T ip43dgoh0e
; real_T hzpktek0ad [ 4 ] ; real_T gmesnt5dc1 ; real_T ljqyvq454t [ 3 ] ;
real_T aejymufvap [ 3 ] ; real_T adoq34mqbb ; real_T iz2k5ya5ta ; real_T
inyi50y04a ; real_T nh5etxksnt ; real_T nuudevirgz ; real_T npg5zigbkh ;
real_T cb55qal4hd ; real_T etztp0mvgy [ 16 ] ; real_T diznsy5jit [ 9 ] ;
int32_T i ; real_T dmy3lnhtsw [ 3 ] ; real_T tmp [ 3 ] ; real_T
b2veiq1ijv_idx_0 ; real_T b2veiq1ijv_idx_1 ; real_T b2veiq1ijv_idx_2 ; real_T
b2veiq1ijv_idx_3 ; pllwawth2v * _rtB ; lfrdjjasqz * _rtP ; isuplyng3d * _rtX
; fw4wgdftov * _rtDW ; _rtDW = ( ( fw4wgdftov * ) ssGetRootDWork ( S ) ) ;
_rtX = ( ( isuplyng3d * ) ssGetContStates ( S ) ) ; _rtP = ( ( lfrdjjasqz * )
ssGetDefaultParam ( S ) ) ; _rtB = ( ( pllwawth2v * ) _ssGetBlockIO ( S ) ) ;
if ( ssIsSampleHit ( S , 1 , tid ) ) { memcpy ( & _rtB -> piujjod4ij [ 0 ] ,
& _rtP -> P_2 [ 0 ] , 9U * sizeof ( real_T ) ) ; } if ( ssIsContinuousTask (
S , tid ) ) { _rtB -> n0qqthcibl [ 0 ] = _rtX -> iidw152zl0 [ 0 ] ; _rtB ->
n0qqthcibl [ 1 ] = _rtX -> iidw152zl0 [ 1 ] ; _rtB -> n0qqthcibl [ 2 ] = _rtX
-> iidw152zl0 [ 2 ] ; for ( i = 0 ; i < 3 ; i ++ ) { fegfldezgg [ i ] = _rtB
-> piujjod4ij [ i + 6 ] * _rtB -> n0qqthcibl [ 2 ] + ( _rtB -> piujjod4ij [ i
+ 3 ] * _rtB -> n0qqthcibl [ 1 ] + _rtB -> piujjod4ij [ i ] * _rtB ->
n0qqthcibl [ 0 ] ) ; } nuudevirgz = _rtP -> P_5 * _rtX -> awllp4kjam ;
npg5zigbkh = _rtP -> P_7 * _rtX -> g11z1syag2 ; cb55qal4hd = _rtP -> P_9 *
_rtX -> bwqacranye ; if ( 0.0 >= _rtP -> P_10 ) { gg2u4if1k1 [ 0 ] = _rtP ->
P_1 * nuudevirgz ; gg2u4if1k1 [ 1 ] = _rtP -> P_1 * npg5zigbkh ; gg2u4if1k1 [
2 ] = _rtP -> P_1 * cb55qal4hd ; } else { gg2u4if1k1 [ 0 ] = _rtP -> P_0 *
nuudevirgz ; gg2u4if1k1 [ 1 ] = _rtP -> P_0 * npg5zigbkh ; gg2u4if1k1 [ 2 ] =
_rtP -> P_0 * cb55qal4hd ; } gg2u4if1k1 [ 0 ] *= _rtP -> P_11 ; gg2u4if1k1 [
1 ] *= _rtP -> P_11 ; gg2u4if1k1 [ 2 ] *= _rtP -> P_11 ; for ( i = 0 ; i < 3
; i ++ ) { _rtB -> fafesurbm0 [ i ] = ( ( _rtB -> piujjod4ij [ i + 3 ] * _rtB
-> n0qqthcibl [ 1 ] + _rtB -> piujjod4ij [ i ] * _rtB -> n0qqthcibl [ 0 ] ) +
_rtB -> piujjod4ij [ i + 6 ] * _rtB -> n0qqthcibl [ 2 ] ) + gg2u4if1k1 [ i ]
; } } if ( ssIsSampleHit ( S , 1 , tid ) ) { for ( i = 0 ; i < 16 ; i ++ ) {
_rtB -> e3xj4w2izs [ i ] = _rtP -> P_12 [ i ] ; _rtB -> gfu3o0binr [ i ] =
_rtP -> P_13 [ i ] ; } } if ( ssIsContinuousTask ( S , tid ) ) { if (
ssGetTaskTime ( S , 0 ) < _rtP -> P_14 ) { h4fomgvj5n = _rtP -> P_15 ; } else
{ h4fomgvj5n = _rtP -> P_16 ; } } if ( ssIsSampleHit ( S , 1 , tid ) ) {
memcpy ( & _rtB -> c0ry2eqmbe [ 0 ] , & _rtP -> P_17 [ 0 ] , sizeof ( real_T
) << 4U ) ; } if ( ssIsContinuousTask ( S , tid ) ) { for ( i = 0 ; i < 16 ;
i ++ ) { if ( h4fomgvj5n >= _rtP -> P_18 ) { etztp0mvgy [ i ] = _rtB ->
gfu3o0binr [ i ] ; } else { etztp0mvgy [ i ] = _rtB -> c0ry2eqmbe [ i ] ; } }
} if ( ssIsSampleHit ( S , 1 , tid ) ) { memcpy ( & _rtB -> desfuzj5lp [ 0 ]
, & _rtP -> P_19 [ 0 ] , sizeof ( real_T ) << 4U ) ; } if (
ssIsContinuousTask ( S , tid ) ) { _rtB -> c3zmymnj3l [ 0 ] = _rtX ->
chvazs3t52 [ 0 ] ; _rtB -> c3zmymnj3l [ 1 ] = _rtX -> chvazs3t52 [ 1 ] ; _rtB
-> c3zmymnj3l [ 2 ] = _rtX -> chvazs3t52 [ 2 ] ; h4fomgvj5n =
muDoubleScalarCos ( _rtB -> c3zmymnj3l [ 0 ] ) ; erv0ccpdoz =
muDoubleScalarCos ( _rtB -> c3zmymnj3l [ 1 ] ) ; jq5l4n1lxp = ssGetT ( S ) ;
jq5l4n1lxp *= _rtP -> P_21 ; _rtB -> lhofsaxnb2 = 0.4 ; _rtB -> ofpcmvzmcw [
0 ] = _rtX -> p4te2bxuor [ 0 ] ; _rtB -> ofpcmvzmcw [ 1 ] = _rtX ->
p4te2bxuor [ 1 ] ; _rtB -> ofpcmvzmcw [ 2 ] = _rtX -> p4te2bxuor [ 2 ] ;
cb55qal4hd = _rtB -> lhofsaxnb2 - _rtB -> ofpcmvzmcw [ 2 ] ; fegfldezgg [ 0 ]
= _rtX -> daauadhpee [ 0 ] ; fegfldezgg [ 1 ] = _rtX -> daauadhpee [ 1 ] ;
fegfldezgg [ 2 ] = _rtX -> daauadhpee [ 2 ] ; _rtB -> h21fxtwche [ 0 ] = _rtX
-> daauadhpee [ 0 ] ; _rtB -> h21fxtwche [ 1 ] = _rtX -> daauadhpee [ 1 ] ;
_rtB -> h21fxtwche [ 2 ] = _rtX -> daauadhpee [ 2 ] ; if ( ( _rtDW ->
pn3mnuc0r2 >= ssGetT ( S ) ) && ( _rtDW -> hjzxgmvdja >= ssGetT ( S ) ) ) {
ip43dgoh0e = 0.0 ; } else { nuudevirgz = _rtDW -> pn3mnuc0r2 ; lastU = &
_rtDW -> kehdx0zpi5 ; if ( _rtDW -> pn3mnuc0r2 < _rtDW -> hjzxgmvdja ) { if (
_rtDW -> hjzxgmvdja < ssGetT ( S ) ) { nuudevirgz = _rtDW -> hjzxgmvdja ;
lastU = & _rtDW -> gqz0n4opfd ; } } else { if ( _rtDW -> pn3mnuc0r2 >= ssGetT
( S ) ) { nuudevirgz = _rtDW -> hjzxgmvdja ; lastU = & _rtDW -> gqz0n4opfd ;
} } ip43dgoh0e = ( _rtB -> lhofsaxnb2 - * lastU ) / ( ssGetT ( S ) -
nuudevirgz ) ; } ip43dgoh0e = ( _rtB -> h21fxtwche [ 2 ] - ip43dgoh0e ) -
_rtP -> P_24 * cb55qal4hd ; _rtB -> exbh02tybg = ( ( ( cb55qal4hd + 9.81 ) -
( 10.0 * cb55qal4hd + ip43dgoh0e ) * 10.0 ) - 0.3 * ip43dgoh0e ) * ( 1.65 /
h4fomgvj5n / erv0ccpdoz ) ; for ( i = 0 ; i < 4 ; i ++ ) { _rtB -> mhxvjc0v0s
[ i ] = 0.0 ; _rtB -> mhxvjc0v0s [ i ] += _rtB -> desfuzj5lp [ i ] * _rtB ->
exbh02tybg ; _rtB -> mhxvjc0v0s [ i ] += _rtB -> desfuzj5lp [ i + 4 ] * _rtB
-> fafesurbm0 [ 0 ] ; _rtB -> mhxvjc0v0s [ i ] += _rtB -> desfuzj5lp [ i + 8
] * _rtB -> fafesurbm0 [ 1 ] ; _rtB -> mhxvjc0v0s [ i ] += _rtB -> desfuzj5lp
[ i + 12 ] * _rtB -> fafesurbm0 [ 2 ] ; } { real_T * * uBuffer = ( real_T * *
) & _rtDW -> ebzxnz4cux . TUbufferPtrs [ 0 ] ; real_T * * tBuffer = ( real_T
* * ) & _rtDW -> ebzxnz4cux . TUbufferPtrs [ 4 ] ; real_T simTime = ssGetT (
S ) ; real_T tMinusDelay ; { int_T i1 ; real_T * y0 = & _rtB -> edskamelsk [
0 ] ; const real_T * u0 = & _rtB -> mhxvjc0v0s [ 0 ] ; int_T * iw_Tail = &
_rtDW -> go5ejd0ffa . Tail [ 0 ] ; int_T * iw_Head = & _rtDW -> go5ejd0ffa .
Head [ 0 ] ; int_T * iw_Last = & _rtDW -> go5ejd0ffa . Last [ 0 ] ; int_T *
iw_CircularBufSize = & _rtDW -> go5ejd0ffa . CircularBufSize [ 0 ] ; for ( i1
= 0 ; i1 < 4 ; i1 ++ ) { tMinusDelay = ( ( _rtP -> P_25 > 0.0 ) ? _rtP ->
P_25 : 0.0 ) ; tMinusDelay = simTime - tMinusDelay ; if ( _rtP -> P_25 == 0.0
) y0 [ i1 ] = u0 [ i1 ] ; else y0 [ i1 ] =
Date0501_QuadModel_acc_rt_TDelayInterpolate ( tMinusDelay , 0.0 , * tBuffer ,
* uBuffer , iw_CircularBufSize [ i1 ] , & iw_Last [ i1 ] , iw_Tail [ i1 ] ,
iw_Head [ i1 ] , _rtP -> P_26 , 0 , ( boolean_T ) ( ssIsMinorTimeStep ( S )
&& ( ssGetTimeOfLastOutput ( S ) == ssGetT ( S ) ) ) ) ; tBuffer ++ ; uBuffer
++ ; } } } ip43dgoh0e = _rtP -> P_27 * _rtB -> edskamelsk [ 0 ] ; cb55qal4hd
= _rtP -> P_28 * _rtB -> edskamelsk [ 1 ] ; erv0ccpdoz = _rtP -> P_29 * _rtB
-> edskamelsk [ 2 ] ; h4fomgvj5n = _rtP -> P_30 * _rtB -> edskamelsk [ 3 ] ;
b2veiq1ijv_idx_0 = ip43dgoh0e ; b2veiq1ijv_idx_1 = cb55qal4hd ;
b2veiq1ijv_idx_2 = erv0ccpdoz ; b2veiq1ijv_idx_3 = h4fomgvj5n ; for ( i = 0 ;
i < 4 ; i ++ ) { npg5zigbkh = etztp0mvgy [ i + 12 ] * h4fomgvj5n + (
etztp0mvgy [ i + 8 ] * erv0ccpdoz + ( etztp0mvgy [ i + 4 ] * cb55qal4hd +
etztp0mvgy [ i ] * ip43dgoh0e ) ) ; hzpktek0ad [ i ] = npg5zigbkh ; } for ( i
= 0 ; i < 4 ; i ++ ) { _rtB -> nxcx0t5nhe [ i ] = 0.0 ; _rtB -> nxcx0t5nhe [
i ] += _rtB -> e3xj4w2izs [ i ] * hzpktek0ad [ 0 ] ; _rtB -> nxcx0t5nhe [ i ]
+= _rtB -> e3xj4w2izs [ i + 4 ] * hzpktek0ad [ 1 ] ; _rtB -> nxcx0t5nhe [ i ]
+= _rtB -> e3xj4w2izs [ i + 8 ] * hzpktek0ad [ 2 ] ; _rtB -> nxcx0t5nhe [ i ]
+= _rtB -> e3xj4w2izs [ i + 12 ] * hzpktek0ad [ 3 ] ; } } if ( ssIsSampleHit
( S , 1 , tid ) ) { ssCallAccelRunBlock ( S , 2 , 41 , SS_CALL_MDL_OUTPUTS )
; ssCallAccelRunBlock ( S , 2 , 42 , SS_CALL_MDL_OUTPUTS ) ; } if (
ssIsContinuousTask ( S , tid ) ) { _rtB -> lufo0jtuu5 = muDoubleScalarCos (
jq5l4n1lxp ) / 4.0 ; ip43dgoh0e = _rtB -> lufo0jtuu5 - _rtB -> ofpcmvzmcw [ 0
] ; if ( ( _rtDW -> e2jci4idik >= ssGetT ( S ) ) && ( _rtDW -> becjzjwo4y >=
ssGetT ( S ) ) ) { cb55qal4hd = 0.0 ; } else { nuudevirgz = _rtDW ->
e2jci4idik ; lastU = & _rtDW -> gtw34zxoo4 ; if ( _rtDW -> e2jci4idik < _rtDW
-> becjzjwo4y ) { if ( _rtDW -> becjzjwo4y < ssGetT ( S ) ) { nuudevirgz =
_rtDW -> becjzjwo4y ; lastU = & _rtDW -> chwddmtxaw ; } } else { if ( _rtDW
-> e2jci4idik >= ssGetT ( S ) ) { nuudevirgz = _rtDW -> becjzjwo4y ; lastU =
& _rtDW -> chwddmtxaw ; } } cb55qal4hd = ( _rtB -> lufo0jtuu5 - * lastU ) / (
ssGetT ( S ) - nuudevirgz ) ; } cb55qal4hd = ( _rtB -> h21fxtwche [ 0 ] -
cb55qal4hd ) - _rtP -> P_31 * ip43dgoh0e ; erv0ccpdoz = ( ( ip43dgoh0e - (
cb55qal4hd + ip43dgoh0e ) ) - 0.1 * cb55qal4hd ) * ( 1.65 / _rtB ->
exbh02tybg ) ; _rtB -> a0g20oizdg = muDoubleScalarSin ( jq5l4n1lxp ) / 4.0 ;
ip43dgoh0e = _rtB -> a0g20oizdg - _rtB -> ofpcmvzmcw [ 1 ] ; if ( ( _rtDW ->
g24m0yoivm >= ssGetT ( S ) ) && ( _rtDW -> o3s5v4luha >= ssGetT ( S ) ) ) {
cb55qal4hd = 0.0 ; } else { nuudevirgz = _rtDW -> g24m0yoivm ; lastU = &
_rtDW -> akzdemflu2 ; if ( _rtDW -> g24m0yoivm < _rtDW -> o3s5v4luha ) { if (
_rtDW -> o3s5v4luha < ssGetT ( S ) ) { nuudevirgz = _rtDW -> o3s5v4luha ;
lastU = & _rtDW -> ewuuifkg0s ; } } else { if ( _rtDW -> g24m0yoivm >= ssGetT
( S ) ) { nuudevirgz = _rtDW -> o3s5v4luha ; lastU = & _rtDW -> ewuuifkg0s ;
} } cb55qal4hd = ( _rtB -> a0g20oizdg - * lastU ) / ( ssGetT ( S ) -
nuudevirgz ) ; } cb55qal4hd = ( _rtB -> h21fxtwche [ 1 ] - cb55qal4hd ) -
_rtP -> P_32 * ip43dgoh0e ; h4fomgvj5n = ( ( ip43dgoh0e - ( cb55qal4hd +
ip43dgoh0e ) ) - 0.1 * cb55qal4hd ) * ( 1.65 / _rtB -> exbh02tybg ) ;
ip43dgoh0e = muDoubleScalarCos ( _rtB -> c3zmymnj3l [ 2 ] ) ; cb55qal4hd =
muDoubleScalarSin ( _rtB -> c3zmymnj3l [ 2 ] ) ; nuudevirgz = cb55qal4hd *
erv0ccpdoz - ip43dgoh0e * h4fomgvj5n ; if ( nuudevirgz > 1.0 ) { nuudevirgz =
1.0 ; } else { if ( nuudevirgz < - 1.0 ) { nuudevirgz = - 1.0 ; } }
gmesnt5dc1 = muDoubleScalarAsin ( nuudevirgz ) ; erv0ccpdoz = ip43dgoh0e *
erv0ccpdoz + cb55qal4hd * h4fomgvj5n ; nuudevirgz = erv0ccpdoz /
muDoubleScalarCos ( gmesnt5dc1 ) ; if ( nuudevirgz > 1.0 ) { nuudevirgz = 1.0
; } else { if ( nuudevirgz < - 1.0 ) { nuudevirgz = - 1.0 ; } } _rtB ->
fjodtin3x5 = muDoubleScalarAsin ( nuudevirgz ) ; jq5l4n1lxp =
muDoubleScalarSin ( jq5l4n1lxp ) ; if ( gmesnt5dc1 > _rtP -> P_33 ) { _rtB ->
hok0j3hzyy [ 0 ] = _rtP -> P_33 ; } else if ( gmesnt5dc1 < _rtP -> P_34 ) {
_rtB -> hok0j3hzyy [ 0 ] = _rtP -> P_34 ; } else { _rtB -> hok0j3hzyy [ 0 ] =
gmesnt5dc1 ; } if ( _rtB -> fjodtin3x5 > _rtP -> P_33 ) { _rtB -> hok0j3hzyy
[ 1 ] = _rtP -> P_33 ; } else if ( _rtB -> fjodtin3x5 < _rtP -> P_34 ) { _rtB
-> hok0j3hzyy [ 1 ] = _rtP -> P_34 ; } else { _rtB -> hok0j3hzyy [ 1 ] = _rtB
-> fjodtin3x5 ; } if ( jq5l4n1lxp > _rtP -> P_33 ) { _rtB -> hok0j3hzyy [ 2 ]
= _rtP -> P_33 ; } else if ( jq5l4n1lxp < _rtP -> P_34 ) { _rtB -> hok0j3hzyy
[ 2 ] = _rtP -> P_34 ; } else { _rtB -> hok0j3hzyy [ 2 ] = jq5l4n1lxp ; } }
if ( ssIsSampleHit ( S , 1 , tid ) ) { ssCallAccelRunBlock ( S , 2 , 65 ,
SS_CALL_MDL_OUTPUTS ) ; ssCallAccelRunBlock ( S , 2 , 66 ,
SS_CALL_MDL_OUTPUTS ) ; ssCallAccelRunBlock ( S , 2 , 67 ,
SS_CALL_MDL_OUTPUTS ) ; ssCallAccelRunBlock ( S , 2 , 68 ,
SS_CALL_MDL_OUTPUTS ) ; ssCallAccelRunBlock ( S , 2 , 69 ,
SS_CALL_MDL_OUTPUTS ) ; ssCallAccelRunBlock ( S , 2 , 70 ,
SS_CALL_MDL_OUTPUTS ) ; } if ( ssIsContinuousTask ( S , tid ) ) { fegfldezgg
[ 0 ] = _rtX -> g1it3hbzj0 [ 0 ] ; fegfldezgg [ 1 ] = _rtX -> g1it3hbzj0 [ 1
] ; fegfldezgg [ 2 ] = _rtX -> g1it3hbzj0 [ 2 ] ; _rtB -> jkjusfctvb [ 0 ] =
_rtX -> g1it3hbzj0 [ 0 ] - _rtB -> n0qqthcibl [ 0 ] ; _rtB -> jkjusfctvb [ 1
] = _rtX -> g1it3hbzj0 [ 1 ] - _rtB -> n0qqthcibl [ 1 ] ; _rtB -> jkjusfctvb
[ 2 ] = _rtX -> g1it3hbzj0 [ 2 ] - _rtB -> n0qqthcibl [ 2 ] ; } if (
ssIsSampleHit ( S , 1 , tid ) ) { ssCallAccelRunBlock ( S , 2 , 74 ,
SS_CALL_MDL_OUTPUTS ) ; } if ( ssIsContinuousTask ( S , tid ) ) { _rtB ->
a4ud5sphga = ( ( _rtP -> P_39 * _rtB -> hok0j3hzyy [ 0 ] - _rtB -> c3zmymnj3l
[ 0 ] ) * _rtP -> P_40 - _rtX -> lypsvf5tuh ) * _rtP -> P_42 ; _rtB ->
eabzjxpk4r = ( ( _rtP -> P_36 * _rtB -> hok0j3hzyy [ 0 ] - _rtB -> c3zmymnj3l
[ 0 ] ) * _rtP -> P_37 + _rtX -> ituablqrnc ) + _rtB -> a4ud5sphga ; _rtB ->
nmsl5lhmkk = ( ( _rtP -> P_46 * _rtB -> hok0j3hzyy [ 1 ] - _rtB -> c3zmymnj3l
[ 1 ] ) * _rtP -> P_47 - _rtX -> gkvrijmbgy ) * _rtP -> P_49 ; _rtB ->
apvzubcyyo = ( ( _rtP -> P_43 * _rtB -> hok0j3hzyy [ 1 ] - _rtB -> c3zmymnj3l
[ 1 ] ) * _rtP -> P_44 + _rtX -> pa1uirjhny ) + _rtB -> nmsl5lhmkk ; _rtB ->
njq2qr1dmi = ( ( _rtP -> P_53 * _rtB -> hok0j3hzyy [ 2 ] - _rtB -> c3zmymnj3l
[ 2 ] ) * _rtP -> P_54 - _rtX -> jgpqgj3dd4 ) * _rtP -> P_56 ; _rtB ->
lub0mt0h4w = ( ( _rtP -> P_50 * _rtB -> hok0j3hzyy [ 2 ] - _rtB -> c3zmymnj3l
[ 2 ] ) * _rtP -> P_51 + _rtX -> jjukormc3p ) + _rtB -> njq2qr1dmi ; } if (
ssIsSampleHit ( S , 1 , tid ) ) { ssCallAccelRunBlock ( S , 2 , 108 ,
SS_CALL_MDL_OUTPUTS ) ; ssCallAccelRunBlock ( S , 2 , 109 ,
SS_CALL_MDL_OUTPUTS ) ; } if ( ssIsSampleHit ( S , 2 , tid ) ) { memcpy ( &
diznsy5jit [ 0 ] , & _rtP -> P_57 [ 0 ] , 9U * sizeof ( real_T ) ) ; } if (
ssIsContinuousTask ( S , tid ) && ssIsSpecialSampleHit ( S , 2 , 0 , tid ) )
{ _rtB -> b54dfnnkoo [ 0 ] = _rtB -> jkjusfctvb [ 0 ] ; _rtB -> b54dfnnkoo [
1 ] = _rtB -> jkjusfctvb [ 1 ] ; _rtB -> b54dfnnkoo [ 2 ] = _rtB ->
jkjusfctvb [ 2 ] ; } if ( ssIsSampleHit ( S , 2 , tid ) ) { for ( i = 0 ; i <
3 ; i ++ ) { dmy3lnhtsw [ i ] = diznsy5jit [ i + 6 ] * _rtB -> b54dfnnkoo [ 2
] + ( diznsy5jit [ i + 3 ] * _rtB -> b54dfnnkoo [ 1 ] + diznsy5jit [ i ] *
_rtB -> b54dfnnkoo [ 0 ] ) ; } _rtB -> atv4umvyhc [ 0 ] = _rtP -> P_58 *
dmy3lnhtsw [ 0 ] ; _rtB -> atv4umvyhc [ 1 ] = _rtP -> P_58 * dmy3lnhtsw [ 1 ]
; _rtB -> atv4umvyhc [ 2 ] = _rtP -> P_58 * dmy3lnhtsw [ 2 ] ; } if (
ssIsSampleHit ( S , 1 , tid ) ) { if ( ssIsSpecialSampleHit ( S , 2 , 1 , tid
) ) { _rtB -> p5bhwepjga [ 0 ] = _rtDW -> nzkmy44stj [ 0 ] ; _rtB ->
p5bhwepjga [ 1 ] = _rtDW -> nzkmy44stj [ 1 ] ; _rtB -> p5bhwepjga [ 2 ] =
_rtDW -> nzkmy44stj [ 2 ] ; } if ( _rtB -> p5bhwepjga [ 0 ] > _rtP -> P_60 )
{ _rtB -> lscax1obsj [ 0 ] = _rtP -> P_60 ; } else if ( _rtB -> p5bhwepjga [
0 ] < _rtP -> P_61 ) { _rtB -> lscax1obsj [ 0 ] = _rtP -> P_61 ; } else {
_rtB -> lscax1obsj [ 0 ] = _rtB -> p5bhwepjga [ 0 ] ; } if ( _rtB ->
p5bhwepjga [ 1 ] > _rtP -> P_60 ) { _rtB -> lscax1obsj [ 1 ] = _rtP -> P_60 ;
} else if ( _rtB -> p5bhwepjga [ 1 ] < _rtP -> P_61 ) { _rtB -> lscax1obsj [
1 ] = _rtP -> P_61 ; } else { _rtB -> lscax1obsj [ 1 ] = _rtB -> p5bhwepjga [
1 ] ; } if ( _rtB -> p5bhwepjga [ 2 ] > _rtP -> P_60 ) { _rtB -> lscax1obsj [
2 ] = _rtP -> P_60 ; } else if ( _rtB -> p5bhwepjga [ 2 ] < _rtP -> P_61 ) {
_rtB -> lscax1obsj [ 2 ] = _rtP -> P_61 ; } else { _rtB -> lscax1obsj [ 2 ] =
_rtB -> p5bhwepjga [ 2 ] ; } memcpy ( & _rtB -> nvfh2rwryr [ 0 ] , & _rtP ->
P_62 [ 0 ] , 9U * sizeof ( real_T ) ) ; } if ( ssIsContinuousTask ( S , tid )
) { ljqyvq454t [ 0 ] = _rtB -> eabzjxpk4r ; ljqyvq454t [ 1 ] = _rtB ->
apvzubcyyo ; ljqyvq454t [ 2 ] = _rtB -> lub0mt0h4w ; for ( i = 0 ; i < 3 ; i
++ ) { aejymufvap [ i ] = _rtB -> nvfh2rwryr [ i + 6 ] * _rtB -> lub0mt0h4w +
( _rtB -> nvfh2rwryr [ i + 3 ] * _rtB -> apvzubcyyo + _rtB -> nvfh2rwryr [ i
] * _rtB -> eabzjxpk4r ) ; } for ( i = 0 ; i < 3 ; i ++ ) { _rtB ->
cnbt5ikdq5 [ i ] = ( _rtB -> lscax1obsj [ i ] - ( ( _rtB -> nvfh2rwryr [ i +
3 ] * _rtB -> apvzubcyyo + _rtB -> nvfh2rwryr [ i ] * _rtB -> eabzjxpk4r ) +
_rtB -> nvfh2rwryr [ i + 6 ] * _rtB -> lub0mt0h4w ) ) + gg2u4if1k1 [ i ] ; }
_rtB -> hgett1hwvh = ( _rtB -> hok0j3hzyy [ 0 ] - _rtB -> c3zmymnj3l [ 0 ] )
* _rtP -> P_63 ; _rtB -> djlmhwhooq = ( _rtB -> hok0j3hzyy [ 1 ] - _rtB ->
c3zmymnj3l [ 1 ] ) * _rtP -> P_64 ; _rtB -> jmv2so5wlt = ( _rtB -> hok0j3hzyy
[ 2 ] - _rtB -> c3zmymnj3l [ 2 ] ) * _rtP -> P_65 ; } if ( ssIsSampleHit ( S
, 1 , tid ) ) { for ( i = 0 ; i < 9 ; i ++ ) { _rtB -> g231mw0l04 [ i ] =
_rtP -> P_66 [ i ] ; _rtB -> cvfhod2fgy [ i ] = _rtP -> P_67 [ i ] ; } } if (
ssIsContinuousTask ( S , tid ) ) { for ( i = 0 ; i < 3 ; i ++ ) { aejymufvap
[ i ] = _rtB -> g231mw0l04 [ i + 6 ] * fegfldezgg [ 2 ] + ( _rtB ->
g231mw0l04 [ i + 3 ] * fegfldezgg [ 1 ] + _rtB -> g231mw0l04 [ i ] *
fegfldezgg [ 0 ] ) ; } gg2u4if1k1 [ 0 ] += _rtB -> lscax1obsj [ 0 ] ;
gg2u4if1k1 [ 1 ] += _rtB -> lscax1obsj [ 1 ] ; npg5zigbkh = gg2u4if1k1 [ 2 ]
+ _rtB -> lscax1obsj [ 2 ] ; for ( i = 0 ; i < 3 ; i ++ ) { ljqyvq454t [ i ]
= _rtB -> cvfhod2fgy [ i + 6 ] * npg5zigbkh + ( _rtB -> cvfhod2fgy [ i + 3 ]
* gg2u4if1k1 [ 1 ] + _rtB -> cvfhod2fgy [ i ] * gg2u4if1k1 [ 0 ] ) ; } for (
i = 0 ; i < 3 ; i ++ ) { tmp [ i ] = _rtB -> cvfhod2fgy [ i + 6 ] *
npg5zigbkh + ( _rtB -> cvfhod2fgy [ i + 3 ] * gg2u4if1k1 [ 1 ] + _rtB ->
cvfhod2fgy [ i ] * gg2u4if1k1 [ 0 ] ) ; } for ( i = 0 ; i < 3 ; i ++ ) {
dmy3lnhtsw [ i ] = _rtB -> g231mw0l04 [ i + 6 ] * fegfldezgg [ 2 ] + ( _rtB
-> g231mw0l04 [ i + 3 ] * fegfldezgg [ 1 ] + _rtB -> g231mw0l04 [ i ] *
fegfldezgg [ 0 ] ) ; } _rtB -> jpqlzrwcim [ 0 ] = tmp [ 0 ] + dmy3lnhtsw [ 0
] ; _rtB -> jpqlzrwcim [ 1 ] = tmp [ 1 ] + dmy3lnhtsw [ 1 ] ; _rtB ->
jpqlzrwcim [ 2 ] = tmp [ 2 ] + dmy3lnhtsw [ 2 ] ; nuudevirgz = erv0ccpdoz /
muDoubleScalarCos ( gmesnt5dc1 ) ; if ( nuudevirgz > 1.0 ) { nuudevirgz = 1.0
; } else { if ( nuudevirgz < - 1.0 ) { nuudevirgz = - 1.0 ; } } _rtB ->
po3tmcjsnj = muDoubleScalarAsin ( nuudevirgz ) ; } if ( ssIsSampleHit ( S , 1
, tid ) ) { ssCallAccelRunBlock ( S , 2 , 133 , SS_CALL_MDL_OUTPUTS ) ; } if
( ssIsContinuousTask ( S , tid ) ) { cb55qal4hd = _rtP -> P_73 * _rtX ->
djg03eirmm ; if ( ssGetTaskTime ( S , 0 ) < _rtP -> P_74 ) { _rtB ->
dwmnl5btul = _rtP -> P_75 ; } else { _rtB -> dwmnl5btul = _rtP -> P_76 ; } if
( ( _rtDW -> nnejqocq0a >= ssGetT ( S ) ) && ( _rtDW -> prputol1xj >= ssGetT
( S ) ) ) { jq5l4n1lxp = 0.0 ; } else { nuudevirgz = _rtDW -> nnejqocq0a ;
lastU = & _rtDW -> nnn0pbu5zk ; if ( _rtDW -> nnejqocq0a < _rtDW ->
prputol1xj ) { if ( _rtDW -> prputol1xj < ssGetT ( S ) ) { nuudevirgz = _rtDW
-> prputol1xj ; lastU = & _rtDW -> itm2s5h21f ; } } else { if ( _rtDW ->
nnejqocq0a >= ssGetT ( S ) ) { nuudevirgz = _rtDW -> prputol1xj ; lastU = &
_rtDW -> itm2s5h21f ; } } jq5l4n1lxp = ( _rtB -> dwmnl5btul - * lastU ) / (
ssGetT ( S ) - nuudevirgz ) ; } if ( jq5l4n1lxp > _rtP -> P_77 ) { jq5l4n1lxp
= _rtP -> P_77 ; } else { if ( jq5l4n1lxp < _rtP -> P_78 ) { jq5l4n1lxp =
_rtP -> P_78 ; } } b2veiq1ijv_idx_0 = ( _rtP -> P_69 * _rtX -> hzmj4zwexc +
_rtB -> nxcx0t5nhe [ 0 ] ) + jq5l4n1lxp ; b2veiq1ijv_idx_1 = ( _rtP -> P_71 *
_rtX -> gdwfbh1eyy + _rtB -> nxcx0t5nhe [ 1 ] ) + 0.0 ; b2veiq1ijv_idx_2 = (
_rtB -> nxcx0t5nhe [ 2 ] + cb55qal4hd ) + jq5l4n1lxp ; b2veiq1ijv_idx_3 = (
_rtB -> nxcx0t5nhe [ 3 ] + cb55qal4hd ) + 0.0 ; } if ( ssIsSampleHit ( S , 1
, tid ) ) { adoq34mqbb = _rtP -> P_82 * _rtDW -> m23rkpyrxj ; iz2k5ya5ta =
_rtP -> P_86 * _rtDW -> kb4sbkwk2c ; inyi50y04a = _rtP -> P_90 * _rtDW ->
hrcvjke4js ; memcpy ( & _rtB -> nmgcamzool [ 0 ] , & _rtP -> P_91 [ 0 ] , 9U
* sizeof ( real_T ) ) ; } if ( ssIsContinuousTask ( S , tid ) ) { for ( i = 0
; i < 3 ; i ++ ) { aejymufvap [ i ] = _rtB -> nmgcamzool [ i + 6 ] * _rtB ->
n0qqthcibl [ 2 ] + ( _rtB -> nmgcamzool [ i + 3 ] * _rtB -> n0qqthcibl [ 1 ]
+ _rtB -> nmgcamzool [ i ] * _rtB -> n0qqthcibl [ 0 ] ) ; } for ( i = 0 ; i <
3 ; i ++ ) { tmp [ i ] = _rtB -> nmgcamzool [ i + 6 ] * _rtB -> n0qqthcibl [
2 ] + ( _rtB -> nmgcamzool [ i + 3 ] * _rtB -> n0qqthcibl [ 1 ] + _rtB ->
nmgcamzool [ i ] * _rtB -> n0qqthcibl [ 0 ] ) ; } jq5l4n1lxp = _rtB ->
n0qqthcibl [ 0 ] * tmp [ 2 ] ; npg5zigbkh = aejymufvap [ 0 ] ; cb55qal4hd =
aejymufvap [ 1 ] ; nuudevirgz = aejymufvap [ 0 ] ; aejymufvap [ 0 ] = _rtB ->
n0qqthcibl [ 1 ] * aejymufvap [ 2 ] - _rtB -> n0qqthcibl [ 2 ] * aejymufvap [
1 ] ; aejymufvap [ 0 ] = _rtP -> P_92 * b2veiq1ijv_idx_1 - aejymufvap [ 0 ] ;
aejymufvap [ 1 ] = _rtP -> P_93 * b2veiq1ijv_idx_2 - ( _rtB -> n0qqthcibl [ 2
] * npg5zigbkh - jq5l4n1lxp ) ; aejymufvap [ 2 ] = b2veiq1ijv_idx_3 - ( _rtB
-> n0qqthcibl [ 0 ] * cb55qal4hd - _rtB -> n0qqthcibl [ 1 ] * nuudevirgz ) ;
} if ( ssIsSampleHit ( S , 1 , tid ) ) { memcpy ( & _rtB -> odclhoxjro [ 0 ]
, & _rtP -> P_94 [ 0 ] , 9U * sizeof ( real_T ) ) ; } if ( ssIsContinuousTask
( S , tid ) ) { for ( i = 0 ; i < 3 ; i ++ ) { ljqyvq454t [ i ] = _rtB ->
odclhoxjro [ i + 6 ] * aejymufvap [ 2 ] + ( _rtB -> odclhoxjro [ i + 3 ] *
aejymufvap [ 1 ] + _rtB -> odclhoxjro [ i ] * aejymufvap [ 0 ] ) ; } } if (
ssIsSampleHit ( S , 1 , tid ) ) { _rtB -> ego0iu3ksj [ 0 ] = _rtP -> P_95 [ 0
] ; _rtB -> ego0iu3ksj [ 1 ] = _rtP -> P_95 [ 1 ] ; _rtB -> ego0iu3ksj [ 2 ]
= _rtP -> P_95 [ 2 ] ; _rtB -> ego0iu3ksj [ 3 ] = _rtP -> P_95 [ 3 ] ; } if (
ssIsContinuousTask ( S , tid ) ) { hzpktek0ad [ 0 ] *= _rtB -> ego0iu3ksj [ 0
] ; hzpktek0ad [ 1 ] *= _rtB -> ego0iu3ksj [ 1 ] ; hzpktek0ad [ 2 ] *= _rtB
-> ego0iu3ksj [ 2 ] ; gmesnt5dc1 = ( ( hzpktek0ad [ 0 ] + hzpktek0ad [ 1 ] )
+ hzpktek0ad [ 2 ] ) + _rtB -> ego0iu3ksj [ 3 ] * hzpktek0ad [ 3 ] ; _rtB ->
fsqy1io23x [ 0 ] = _rtP -> P_96 * gmesnt5dc1 * _rtB -> n0qqthcibl [ 1 ] +
ljqyvq454t [ 0 ] ; _rtB -> fsqy1io23x [ 1 ] = _rtP -> P_97 * gmesnt5dc1 *
_rtB -> n0qqthcibl [ 0 ] + ljqyvq454t [ 1 ] ; _rtB -> fsqy1io23x [ 2 ] =
ljqyvq454t [ 2 ] + 0.0 ; } if ( ssIsSampleHit ( S , 1 , tid ) ) { _rtB ->
jycqszl3zk = _rtP -> P_98 * iz2k5ya5ta ; _rtB -> pc0qkleume = _rtP -> P_99 *
adoq34mqbb ; _rtB -> l5y5ceeztf = _rtP -> P_100 * inyi50y04a ; } if (
ssIsContinuousTask ( S , tid ) ) { gmesnt5dc1 = muDoubleScalarCos ( _rtB ->
c3zmymnj3l [ 0 ] ) ; ip43dgoh0e = muDoubleScalarSin ( _rtB -> c3zmymnj3l [ 1
] ) ; cb55qal4hd = muDoubleScalarCos ( _rtB -> c3zmymnj3l [ 2 ] ) ;
jq5l4n1lxp = muDoubleScalarSin ( _rtB -> c3zmymnj3l [ 0 ] ) ; erv0ccpdoz =
muDoubleScalarSin ( _rtB -> c3zmymnj3l [ 2 ] ) ; _rtB -> llxismvqlq = (
gmesnt5dc1 * ip43dgoh0e * cb55qal4hd + jq5l4n1lxp * erv0ccpdoz ) *
b2veiq1ijv_idx_0 * _rtP -> P_101 ; _rtB -> g5w3bgcnkd = ( gmesnt5dc1 *
ip43dgoh0e * erv0ccpdoz - jq5l4n1lxp * cb55qal4hd ) * b2veiq1ijv_idx_0 * _rtP
-> P_102 ; nh5etxksnt = gmesnt5dc1 * muDoubleScalarCos ( _rtB -> c3zmymnj3l [
1 ] ) * b2veiq1ijv_idx_0 * _rtP -> P_103 ; } if ( ssIsSampleHit ( S , 1 , tid
) ) { _rtB -> oinnib4j50 = _rtP -> P_104 ; } if ( ssIsContinuousTask ( S ,
tid ) ) { _rtB -> bavopr1xos = _rtB -> oinnib4j50 + nh5etxksnt ; } if (
ssIsSampleHit ( S , 1 , tid ) ) { ssCallAccelRunBlock ( S , 2 , 193 ,
SS_CALL_MDL_OUTPUTS ) ; ssCallAccelRunBlock ( S , 2 , 194 ,
SS_CALL_MDL_OUTPUTS ) ; } UNUSED_PARAMETER ( tid ) ; }
#define MDL_UPDATE
static void mdlUpdate ( SimStruct * S , int_T tid ) { real_T * lastU ; char_T
* sErr ; pllwawth2v * _rtB ; lfrdjjasqz * _rtP ; fw4wgdftov * _rtDW ; _rtDW =
( ( fw4wgdftov * ) ssGetRootDWork ( S ) ) ; _rtP = ( ( lfrdjjasqz * )
ssGetDefaultParam ( S ) ) ; _rtB = ( ( pllwawth2v * ) _ssGetBlockIO ( S ) ) ;
if ( ssIsContinuousTask ( S , tid ) ) { if ( _rtDW -> pn3mnuc0r2 == ( rtInf )
) { _rtDW -> pn3mnuc0r2 = ssGetT ( S ) ; lastU = & _rtDW -> kehdx0zpi5 ; }
else if ( _rtDW -> hjzxgmvdja == ( rtInf ) ) { _rtDW -> hjzxgmvdja = ssGetT (
S ) ; lastU = & _rtDW -> gqz0n4opfd ; } else if ( _rtDW -> pn3mnuc0r2 < _rtDW
-> hjzxgmvdja ) { _rtDW -> pn3mnuc0r2 = ssGetT ( S ) ; lastU = & _rtDW ->
kehdx0zpi5 ; } else { _rtDW -> hjzxgmvdja = ssGetT ( S ) ; lastU = & _rtDW ->
gqz0n4opfd ; } * lastU = _rtB -> lhofsaxnb2 ; { real_T * * uBuffer = ( real_T
* * ) & _rtDW -> ebzxnz4cux . TUbufferPtrs [ 0 ] ; real_T * * tBuffer = (
real_T * * ) & _rtDW -> ebzxnz4cux . TUbufferPtrs [ 4 ] ; real_T simTime =
ssGetT ( S ) ; _rtDW -> go5ejd0ffa . Head [ 0 ] = ( ( _rtDW -> go5ejd0ffa .
Head [ 0 ] < ( _rtDW -> go5ejd0ffa . CircularBufSize [ 0 ] - 1 ) ) ? ( _rtDW
-> go5ejd0ffa . Head [ 0 ] + 1 ) : 0 ) ; if ( _rtDW -> go5ejd0ffa . Head [ 0
] == _rtDW -> go5ejd0ffa . Tail [ 0 ] ) { if ( !
Date0501_QuadModel_acc_rt_TDelayUpdateTailOrGrowBuf ( & _rtDW -> go5ejd0ffa .
CircularBufSize [ 0 ] , & _rtDW -> go5ejd0ffa . Tail [ 0 ] , & _rtDW ->
go5ejd0ffa . Head [ 0 ] , & _rtDW -> go5ejd0ffa . Last [ 0 ] , simTime - _rtP
-> P_25 , tBuffer , uBuffer , ( NULL ) , ( boolean_T ) 0 , false , & _rtDW ->
go5ejd0ffa . MaxNewBufSize ) ) { ssSetErrorStatus ( S ,
"tdelay memory allocation error" ) ; return ; } } ( * tBuffer ++ ) [ _rtDW ->
go5ejd0ffa . Head [ 0 ] ] = simTime ; ( * uBuffer ++ ) [ _rtDW -> go5ejd0ffa
. Head [ 0 ] ] = _rtB -> mhxvjc0v0s [ 0 ] ; _rtDW -> go5ejd0ffa . Head [ 1 ]
= ( ( _rtDW -> go5ejd0ffa . Head [ 1 ] < ( _rtDW -> go5ejd0ffa .
CircularBufSize [ 1 ] - 1 ) ) ? ( _rtDW -> go5ejd0ffa . Head [ 1 ] + 1 ) : 0
) ; if ( _rtDW -> go5ejd0ffa . Head [ 1 ] == _rtDW -> go5ejd0ffa . Tail [ 1 ]
) { if ( ! Date0501_QuadModel_acc_rt_TDelayUpdateTailOrGrowBuf ( & _rtDW ->
go5ejd0ffa . CircularBufSize [ 1 ] , & _rtDW -> go5ejd0ffa . Tail [ 1 ] , &
_rtDW -> go5ejd0ffa . Head [ 1 ] , & _rtDW -> go5ejd0ffa . Last [ 1 ] ,
simTime - _rtP -> P_25 , tBuffer , uBuffer , ( NULL ) , ( boolean_T ) 0 ,
false , & _rtDW -> go5ejd0ffa . MaxNewBufSize ) ) { ssSetErrorStatus ( S ,
"tdelay memory allocation error" ) ; return ; } } ( * tBuffer ++ ) [ _rtDW ->
go5ejd0ffa . Head [ 1 ] ] = simTime ; ( * uBuffer ++ ) [ _rtDW -> go5ejd0ffa
. Head [ 1 ] ] = _rtB -> mhxvjc0v0s [ 1 ] ; _rtDW -> go5ejd0ffa . Head [ 2 ]
= ( ( _rtDW -> go5ejd0ffa . Head [ 2 ] < ( _rtDW -> go5ejd0ffa .
CircularBufSize [ 2 ] - 1 ) ) ? ( _rtDW -> go5ejd0ffa . Head [ 2 ] + 1 ) : 0
) ; if ( _rtDW -> go5ejd0ffa . Head [ 2 ] == _rtDW -> go5ejd0ffa . Tail [ 2 ]
) { if ( ! Date0501_QuadModel_acc_rt_TDelayUpdateTailOrGrowBuf ( & _rtDW ->
go5ejd0ffa . CircularBufSize [ 2 ] , & _rtDW -> go5ejd0ffa . Tail [ 2 ] , &
_rtDW -> go5ejd0ffa . Head [ 2 ] , & _rtDW -> go5ejd0ffa . Last [ 2 ] ,
simTime - _rtP -> P_25 , tBuffer , uBuffer , ( NULL ) , ( boolean_T ) 0 ,
false , & _rtDW -> go5ejd0ffa . MaxNewBufSize ) ) { ssSetErrorStatus ( S ,
"tdelay memory allocation error" ) ; return ; } } ( * tBuffer ++ ) [ _rtDW ->
go5ejd0ffa . Head [ 2 ] ] = simTime ; ( * uBuffer ++ ) [ _rtDW -> go5ejd0ffa
. Head [ 2 ] ] = _rtB -> mhxvjc0v0s [ 2 ] ; _rtDW -> go5ejd0ffa . Head [ 3 ]
= ( ( _rtDW -> go5ejd0ffa . Head [ 3 ] < ( _rtDW -> go5ejd0ffa .
CircularBufSize [ 3 ] - 1 ) ) ? ( _rtDW -> go5ejd0ffa . Head [ 3 ] + 1 ) : 0
) ; if ( _rtDW -> go5ejd0ffa . Head [ 3 ] == _rtDW -> go5ejd0ffa . Tail [ 3 ]
) { if ( ! Date0501_QuadModel_acc_rt_TDelayUpdateTailOrGrowBuf ( & _rtDW ->
go5ejd0ffa . CircularBufSize [ 3 ] , & _rtDW -> go5ejd0ffa . Tail [ 3 ] , &
_rtDW -> go5ejd0ffa . Head [ 3 ] , & _rtDW -> go5ejd0ffa . Last [ 3 ] ,
simTime - _rtP -> P_25 , tBuffer , uBuffer , ( NULL ) , ( boolean_T ) 0 ,
false , & _rtDW -> go5ejd0ffa . MaxNewBufSize ) ) { ssSetErrorStatus ( S ,
"tdelay memory allocation error" ) ; return ; } } ( * tBuffer ) [ _rtDW ->
go5ejd0ffa . Head [ 3 ] ] = simTime ; ( * uBuffer ) [ _rtDW -> go5ejd0ffa .
Head [ 3 ] ] = _rtB -> mhxvjc0v0s [ 3 ] ; } } if ( ssIsContinuousTask ( S ,
tid ) ) { if ( _rtDW -> e2jci4idik == ( rtInf ) ) { _rtDW -> e2jci4idik =
ssGetT ( S ) ; lastU = & _rtDW -> gtw34zxoo4 ; } else if ( _rtDW ->
becjzjwo4y == ( rtInf ) ) { _rtDW -> becjzjwo4y = ssGetT ( S ) ; lastU = &
_rtDW -> chwddmtxaw ; } else if ( _rtDW -> e2jci4idik < _rtDW -> becjzjwo4y )
{ _rtDW -> e2jci4idik = ssGetT ( S ) ; lastU = & _rtDW -> gtw34zxoo4 ; } else
{ _rtDW -> becjzjwo4y = ssGetT ( S ) ; lastU = & _rtDW -> chwddmtxaw ; } *
lastU = _rtB -> lufo0jtuu5 ; if ( _rtDW -> g24m0yoivm == ( rtInf ) ) { _rtDW
-> g24m0yoivm = ssGetT ( S ) ; lastU = & _rtDW -> akzdemflu2 ; } else if (
_rtDW -> o3s5v4luha == ( rtInf ) ) { _rtDW -> o3s5v4luha = ssGetT ( S ) ;
lastU = & _rtDW -> ewuuifkg0s ; } else if ( _rtDW -> g24m0yoivm < _rtDW ->
o3s5v4luha ) { _rtDW -> g24m0yoivm = ssGetT ( S ) ; lastU = & _rtDW ->
akzdemflu2 ; } else { _rtDW -> o3s5v4luha = ssGetT ( S ) ; lastU = & _rtDW ->
ewuuifkg0s ; } * lastU = _rtB -> a0g20oizdg ; } if ( ssIsSampleHit ( S , 1 ,
tid ) ) { sErr = GetErrorBuffer ( & _rtDW -> byjyo5m5d4 [ 0U ] ) ;
LibUpdate_Network ( & _rtDW -> byjyo5m5d4 [ 0U ] , & _rtB -> pilq1cnd5o [ 0U
] , 48 ) ; if ( * sErr != 0 ) { ssSetErrorStatus ( S , sErr ) ;
ssSetStopRequested ( S , 1 ) ; } } if ( ssIsSampleHit ( S , 2 , tid ) ) {
_rtDW -> nzkmy44stj [ 0 ] = _rtB -> atv4umvyhc [ 0 ] ; _rtDW -> nzkmy44stj [
1 ] = _rtB -> atv4umvyhc [ 1 ] ; _rtDW -> nzkmy44stj [ 2 ] = _rtB ->
atv4umvyhc [ 2 ] ; } if ( ssIsContinuousTask ( S , tid ) ) { if ( _rtDW ->
nnejqocq0a == ( rtInf ) ) { _rtDW -> nnejqocq0a = ssGetT ( S ) ; lastU = &
_rtDW -> nnn0pbu5zk ; } else if ( _rtDW -> prputol1xj == ( rtInf ) ) { _rtDW
-> prputol1xj = ssGetT ( S ) ; lastU = & _rtDW -> itm2s5h21f ; } else if (
_rtDW -> nnejqocq0a < _rtDW -> prputol1xj ) { _rtDW -> nnejqocq0a = ssGetT (
S ) ; lastU = & _rtDW -> nnn0pbu5zk ; } else { _rtDW -> prputol1xj = ssGetT (
S ) ; lastU = & _rtDW -> itm2s5h21f ; } * lastU = _rtB -> dwmnl5btul ; } if (
ssIsSampleHit ( S , 1 , tid ) ) { _rtDW -> m23rkpyrxj =
rt_nrand_Upu32_Yd_f_pw_snf ( & _rtDW -> hugx1ib0cr ) * _rtP -> P_80 + _rtP ->
P_79 ; _rtDW -> kb4sbkwk2c = rt_nrand_Upu32_Yd_f_pw_snf ( & _rtDW ->
nqrt1fyklg ) * _rtP -> P_84 + _rtP -> P_83 ; _rtDW -> hrcvjke4js =
rt_nrand_Upu32_Yd_f_pw_snf ( & _rtDW -> btnlgyxfg3 ) * _rtP -> P_88 + _rtP ->
P_87 ; } UNUSED_PARAMETER ( tid ) ; }
#define MDL_DERIVATIVES
static void mdlDerivatives ( SimStruct * S ) { pllwawth2v * _rtB ; lfrdjjasqz
* _rtP ; isuplyng3d * _rtX ; ns1zbfbemf * _rtXdot ; _rtXdot = ( ( ns1zbfbemf
* ) ssGetdX ( S ) ) ; _rtX = ( ( isuplyng3d * ) ssGetContStates ( S ) ) ;
_rtP = ( ( lfrdjjasqz * ) ssGetDefaultParam ( S ) ) ; _rtB = ( ( pllwawth2v *
) _ssGetBlockIO ( S ) ) ; _rtXdot -> iidw152zl0 [ 0 ] = _rtB -> fsqy1io23x [
0 ] ; _rtXdot -> iidw152zl0 [ 1 ] = _rtB -> fsqy1io23x [ 1 ] ; _rtXdot ->
iidw152zl0 [ 2 ] = _rtB -> fsqy1io23x [ 2 ] ; _rtXdot -> awllp4kjam = 0.0 ;
_rtXdot -> awllp4kjam += _rtP -> P_4 * _rtX -> awllp4kjam ; _rtXdot ->
awllp4kjam += _rtB -> cnbt5ikdq5 [ 0 ] ; _rtXdot -> g11z1syag2 = 0.0 ;
_rtXdot -> g11z1syag2 += _rtP -> P_6 * _rtX -> g11z1syag2 ; _rtXdot ->
g11z1syag2 += _rtB -> cnbt5ikdq5 [ 1 ] ; _rtXdot -> bwqacranye = 0.0 ;
_rtXdot -> bwqacranye += _rtP -> P_8 * _rtX -> bwqacranye ; _rtXdot ->
bwqacranye += _rtB -> cnbt5ikdq5 [ 2 ] ; _rtXdot -> chvazs3t52 [ 0 ] = _rtB
-> n0qqthcibl [ 0 ] ; _rtXdot -> chvazs3t52 [ 1 ] = _rtB -> n0qqthcibl [ 1 ]
; _rtXdot -> chvazs3t52 [ 2 ] = _rtB -> n0qqthcibl [ 2 ] ; _rtXdot ->
p4te2bxuor [ 0 ] = _rtB -> h21fxtwche [ 0 ] ; _rtXdot -> p4te2bxuor [ 1 ] =
_rtB -> h21fxtwche [ 1 ] ; _rtXdot -> p4te2bxuor [ 2 ] = _rtB -> h21fxtwche [
2 ] ; _rtXdot -> daauadhpee [ 0 ] = _rtB -> llxismvqlq ; _rtXdot ->
daauadhpee [ 1 ] = _rtB -> g5w3bgcnkd ; _rtXdot -> daauadhpee [ 2 ] = _rtB ->
bavopr1xos ; _rtXdot -> g1it3hbzj0 [ 0 ] = _rtB -> jpqlzrwcim [ 0 ] ; _rtXdot
-> g1it3hbzj0 [ 1 ] = _rtB -> jpqlzrwcim [ 1 ] ; _rtXdot -> g1it3hbzj0 [ 2 ]
= _rtB -> jpqlzrwcim [ 2 ] ; _rtXdot -> ituablqrnc = _rtB -> hgett1hwvh ;
_rtXdot -> lypsvf5tuh = _rtB -> a4ud5sphga ; _rtXdot -> pa1uirjhny = _rtB ->
djlmhwhooq ; _rtXdot -> gkvrijmbgy = _rtB -> nmsl5lhmkk ; _rtXdot ->
jjukormc3p = _rtB -> jmv2so5wlt ; _rtXdot -> jgpqgj3dd4 = _rtB -> njq2qr1dmi
; _rtXdot -> hzmj4zwexc = 0.0 ; _rtXdot -> hzmj4zwexc += _rtP -> P_68 * _rtX
-> hzmj4zwexc ; _rtXdot -> hzmj4zwexc += _rtB -> pc0qkleume ; _rtXdot ->
gdwfbh1eyy = 0.0 ; _rtXdot -> gdwfbh1eyy += _rtP -> P_70 * _rtX -> gdwfbh1eyy
; _rtXdot -> gdwfbh1eyy += _rtB -> jycqszl3zk ; _rtXdot -> djg03eirmm = 0.0 ;
_rtXdot -> djg03eirmm += _rtP -> P_72 * _rtX -> djg03eirmm ; _rtXdot ->
djg03eirmm += _rtB -> l5y5ceeztf ; } static void mdlInitializeSizes (
SimStruct * S ) { ssSetChecksumVal ( S , 0 , 1803018131U ) ; ssSetChecksumVal
( S , 1 , 3199367205U ) ; ssSetChecksumVal ( S , 2 , 1798073909U ) ;
ssSetChecksumVal ( S , 3 , 284109315U ) ; { mxArray * slVerStructMat = NULL ;
mxArray * slStrMat = mxCreateString ( "simulink" ) ; char slVerChar [ 10 ] ;
int status = mexCallMATLAB ( 1 , & slVerStructMat , 1 , & slStrMat , "ver" )
; if ( status == 0 ) { mxArray * slVerMat = mxGetField ( slVerStructMat , 0 ,
"Version" ) ; if ( slVerMat == NULL ) { status = 1 ; } else { status =
mxGetString ( slVerMat , slVerChar , 10 ) ; } } mxDestroyArray ( slStrMat ) ;
mxDestroyArray ( slVerStructMat ) ; if ( ( status == 1 ) || ( strcmp (
slVerChar , "8.5" ) != 0 ) ) { return ; } } ssSetOptions ( S ,
SS_OPTION_EXCEPTION_FREE_CODE ) ; if ( ssGetSizeofDWork ( S ) != sizeof (
fw4wgdftov ) ) { ssSetErrorStatus ( S ,
"Unexpected error: Internal DWork sizes do "
"not match for accelerator mex file." ) ; } if ( ssGetSizeofGlobalBlockIO ( S
) != sizeof ( pllwawth2v ) ) { ssSetErrorStatus ( S ,
"Unexpected error: Internal BlockIO sizes do "
"not match for accelerator mex file." ) ; } { int ssSizeofParams ;
ssGetSizeofParams ( S , & ssSizeofParams ) ; if ( ssSizeofParams != sizeof (
lfrdjjasqz ) ) { static char msg [ 256 ] ; sprintf ( msg ,
"Unexpected error: Internal Parameters sizes do "
"not match for accelerator mex file." ) ; } } _ssSetDefaultParam ( S , (
real_T * ) & lp5yphtiof ) ; rt_InitInfAndNaN ( sizeof ( real_T ) ) ; } static
void mdlInitializeSampleTimes ( SimStruct * S ) { } static void mdlTerminate
( SimStruct * S ) { }
#include "simulink.c"
