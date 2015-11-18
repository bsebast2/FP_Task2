#include "__cf_Date0501_QuadModel_Thiago.h"
#include <math.h>
#include "Date0501_QuadModel_Thiago_acc.h"
#include "Date0501_QuadModel_Thiago_acc_private.h"
#include <stdio.h>
#include "simstruc.h"
#include "fixedpoint.h"
#define CodeFormat S-Function
#define AccDefine1 Accelerator_S-Function
#ifndef __RTW_UTFREE__  
extern void * utMalloc ( size_t ) ; extern void utFree ( void * ) ;
#endif
boolean_T Date0501_QuadModel_Thiago_acc_rt_TDelayUpdateTailOrGrowBuf ( int_T
* bufSzPtr , int_T * tailPtr , int_T * headPtr , int_T * lastPtr , real_T
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
return ( true ) ; } real_T Date0501_QuadModel_Thiago_acc_rt_TDelayInterpolate
( real_T tMinusDelay , real_T tStart , real_T * tBuf , real_T * uBuf , int_T
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
real_T gshcvxpyxd ; real_T * lastU ; real_T fegfldezgg [ 3 ] ; real_T
gg2u4if1k1 [ 3 ] ; real_T h4fomgvj5n ; real_T erv0ccpdoz ; real_T ip43dgoh0e
; real_T hzpktek0ad [ 4 ] ; real_T gmesnt5dc1 ; real_T ljqyvq454t [ 3 ] ;
real_T aejymufvap [ 3 ] ; real_T adoq34mqbb ; real_T iz2k5ya5ta ; real_T
inyi50y04a ; real_T nh5etxksnt ; real_T nuudevirgz ; real_T npg5zigbkh ;
real_T cb55qal4hd ; real_T etztp0mvgy [ 16 ] ; real_T diznsy5jit [ 9 ] ;
int32_T i ; real_T fvvgxc4onr [ 3 ] ; real_T tmp [ 3 ] ; real_T
a4qrcc3br4_idx_0 ; real_T a4qrcc3br4_idx_1 ; real_T a4qrcc3br4_idx_2 ; real_T
a4qrcc3br4_idx_3 ; gsegsuzktn * _rtB ; pcuunorfmd * _rtP ; pvlbxo1noz * _rtX
; pfm4mif5gq * _rtDW ; _rtDW = ( ( pfm4mif5gq * ) ssGetRootDWork ( S ) ) ;
_rtX = ( ( pvlbxo1noz * ) ssGetContStates ( S ) ) ; _rtP = ( ( pcuunorfmd * )
ssGetDefaultParam ( S ) ) ; _rtB = ( ( gsegsuzktn * ) _ssGetBlockIO ( S ) ) ;
if ( ssIsSampleHit ( S , 1 , tid ) ) { memcpy ( & _rtB -> inlmvodqay [ 0 ] ,
& _rtP -> P_2 [ 0 ] , 9U * sizeof ( real_T ) ) ; } if ( ssIsContinuousTask (
S , tid ) ) { _rtB -> maivwxyydn [ 0 ] = _rtX -> mdqhadbmzi [ 0 ] ; _rtB ->
maivwxyydn [ 1 ] = _rtX -> mdqhadbmzi [ 1 ] ; _rtB -> maivwxyydn [ 2 ] = _rtX
-> mdqhadbmzi [ 2 ] ; for ( i = 0 ; i < 3 ; i ++ ) { fegfldezgg [ i ] = _rtB
-> inlmvodqay [ i + 6 ] * _rtB -> maivwxyydn [ 2 ] + ( _rtB -> inlmvodqay [ i
+ 3 ] * _rtB -> maivwxyydn [ 1 ] + _rtB -> inlmvodqay [ i ] * _rtB ->
maivwxyydn [ 0 ] ) ; } nuudevirgz = _rtP -> P_5 * _rtX -> kswdq12pcl ;
npg5zigbkh = _rtP -> P_7 * _rtX -> bpy1pus01d ; cb55qal4hd = _rtP -> P_9 *
_rtX -> pqy1xvwazu ; if ( 0.0 >= _rtP -> P_10 ) { gg2u4if1k1 [ 0 ] = _rtP ->
P_1 * nuudevirgz ; gg2u4if1k1 [ 1 ] = _rtP -> P_1 * npg5zigbkh ; gg2u4if1k1 [
2 ] = _rtP -> P_1 * cb55qal4hd ; } else { gg2u4if1k1 [ 0 ] = _rtP -> P_0 *
nuudevirgz ; gg2u4if1k1 [ 1 ] = _rtP -> P_0 * npg5zigbkh ; gg2u4if1k1 [ 2 ] =
_rtP -> P_0 * cb55qal4hd ; } gg2u4if1k1 [ 0 ] *= _rtP -> P_11 ; gg2u4if1k1 [
1 ] *= _rtP -> P_11 ; gg2u4if1k1 [ 2 ] *= _rtP -> P_11 ; for ( i = 0 ; i < 3
; i ++ ) { _rtB -> d0yxlwjjin [ i ] = ( ( _rtB -> inlmvodqay [ i + 3 ] * _rtB
-> maivwxyydn [ 1 ] + _rtB -> inlmvodqay [ i ] * _rtB -> maivwxyydn [ 0 ] ) +
_rtB -> inlmvodqay [ i + 6 ] * _rtB -> maivwxyydn [ 2 ] ) + gg2u4if1k1 [ i ]
; } } if ( ssIsSampleHit ( S , 1 , tid ) ) { for ( i = 0 ; i < 16 ; i ++ ) {
_rtB -> i51wezcj32 [ i ] = _rtP -> P_12 [ i ] ; _rtB -> joxpnqy42q [ i ] =
_rtP -> P_13 [ i ] ; } } if ( ssIsContinuousTask ( S , tid ) ) { if (
ssGetTaskTime ( S , 0 ) < _rtP -> P_14 ) { h4fomgvj5n = _rtP -> P_15 ; } else
{ h4fomgvj5n = _rtP -> P_16 ; } } if ( ssIsSampleHit ( S , 1 , tid ) ) {
memcpy ( & _rtB -> p5xiqvzfqn [ 0 ] , & _rtP -> P_17 [ 0 ] , sizeof ( real_T
) << 4U ) ; } if ( ssIsContinuousTask ( S , tid ) ) { for ( i = 0 ; i < 16 ;
i ++ ) { if ( h4fomgvj5n >= _rtP -> P_18 ) { etztp0mvgy [ i ] = _rtB ->
joxpnqy42q [ i ] ; } else { etztp0mvgy [ i ] = _rtB -> p5xiqvzfqn [ i ] ; } }
} if ( ssIsSampleHit ( S , 1 , tid ) ) { memcpy ( & _rtB -> dkvdt2n2h4 [ 0 ]
, & _rtP -> P_19 [ 0 ] , sizeof ( real_T ) << 4U ) ; } if (
ssIsContinuousTask ( S , tid ) ) { _rtB -> fl4bn1wmqu [ 0 ] = _rtX ->
dwor1dwwuw [ 0 ] ; _rtB -> fl4bn1wmqu [ 1 ] = _rtX -> dwor1dwwuw [ 1 ] ; _rtB
-> fl4bn1wmqu [ 2 ] = _rtX -> dwor1dwwuw [ 2 ] ; h4fomgvj5n =
muDoubleScalarCos ( _rtB -> fl4bn1wmqu [ 0 ] ) ; erv0ccpdoz =
muDoubleScalarCos ( _rtB -> fl4bn1wmqu [ 1 ] ) ; gshcvxpyxd = ssGetT ( S ) ;
gshcvxpyxd *= _rtP -> P_21 ; _rtB -> b3px43da1a = 0.4 ; _rtB -> crqrmutysw [
0 ] = _rtX -> mseabyslip [ 0 ] ; _rtB -> crqrmutysw [ 1 ] = _rtX ->
mseabyslip [ 1 ] ; _rtB -> crqrmutysw [ 2 ] = _rtX -> mseabyslip [ 2 ] ;
cb55qal4hd = _rtB -> b3px43da1a - _rtB -> crqrmutysw [ 2 ] ; fegfldezgg [ 0 ]
= _rtX -> mfpumpgq5t [ 0 ] ; fegfldezgg [ 1 ] = _rtX -> mfpumpgq5t [ 1 ] ;
fegfldezgg [ 2 ] = _rtX -> mfpumpgq5t [ 2 ] ; _rtB -> np4c234dzr [ 0 ] = _rtX
-> mfpumpgq5t [ 0 ] ; _rtB -> np4c234dzr [ 1 ] = _rtX -> mfpumpgq5t [ 1 ] ;
_rtB -> np4c234dzr [ 2 ] = _rtX -> mfpumpgq5t [ 2 ] ; if ( ( _rtDW ->
kmpk31urxn >= ssGetT ( S ) ) && ( _rtDW -> ogdw1xlqmp >= ssGetT ( S ) ) ) {
ip43dgoh0e = 0.0 ; } else { nuudevirgz = _rtDW -> kmpk31urxn ; lastU = &
_rtDW -> hphjcnocig ; if ( _rtDW -> kmpk31urxn < _rtDW -> ogdw1xlqmp ) { if (
_rtDW -> ogdw1xlqmp < ssGetT ( S ) ) { nuudevirgz = _rtDW -> ogdw1xlqmp ;
lastU = & _rtDW -> jjecbj5saz ; } } else { if ( _rtDW -> kmpk31urxn >= ssGetT
( S ) ) { nuudevirgz = _rtDW -> ogdw1xlqmp ; lastU = & _rtDW -> jjecbj5saz ;
} } ip43dgoh0e = ( _rtB -> b3px43da1a - * lastU ) / ( ssGetT ( S ) -
nuudevirgz ) ; } ip43dgoh0e = ( _rtB -> np4c234dzr [ 2 ] - ip43dgoh0e ) -
_rtP -> P_24 * cb55qal4hd ; _rtB -> c4oqtnupxw = ( ( ( cb55qal4hd + 9.81 ) -
( 10.0 * cb55qal4hd + ip43dgoh0e ) * 10.0 ) - 0.3 * ip43dgoh0e ) * ( 1.65 /
h4fomgvj5n / erv0ccpdoz ) ; for ( i = 0 ; i < 4 ; i ++ ) { _rtB -> js2p3cmsfo
[ i ] = 0.0 ; _rtB -> js2p3cmsfo [ i ] += _rtB -> dkvdt2n2h4 [ i ] * _rtB ->
c4oqtnupxw ; _rtB -> js2p3cmsfo [ i ] += _rtB -> dkvdt2n2h4 [ i + 4 ] * _rtB
-> d0yxlwjjin [ 0 ] ; _rtB -> js2p3cmsfo [ i ] += _rtB -> dkvdt2n2h4 [ i + 8
] * _rtB -> d0yxlwjjin [ 1 ] ; _rtB -> js2p3cmsfo [ i ] += _rtB -> dkvdt2n2h4
[ i + 12 ] * _rtB -> d0yxlwjjin [ 2 ] ; } { real_T * * uBuffer = ( real_T * *
) & _rtDW -> nulonpacio . TUbufferPtrs [ 0 ] ; real_T * * tBuffer = ( real_T
* * ) & _rtDW -> nulonpacio . TUbufferPtrs [ 4 ] ; real_T simTime = ssGetT (
S ) ; real_T tMinusDelay ; { int_T i1 ; real_T * y0 = & _rtB -> az2wtq5aqx [
0 ] ; const real_T * u0 = & _rtB -> js2p3cmsfo [ 0 ] ; int_T * iw_Tail = &
_rtDW -> ijfvsgcbuh . Tail [ 0 ] ; int_T * iw_Head = & _rtDW -> ijfvsgcbuh .
Head [ 0 ] ; int_T * iw_Last = & _rtDW -> ijfvsgcbuh . Last [ 0 ] ; int_T *
iw_CircularBufSize = & _rtDW -> ijfvsgcbuh . CircularBufSize [ 0 ] ; for ( i1
= 0 ; i1 < 4 ; i1 ++ ) { tMinusDelay = ( ( _rtP -> P_25 > 0.0 ) ? _rtP ->
P_25 : 0.0 ) ; tMinusDelay = simTime - tMinusDelay ; if ( _rtP -> P_25 == 0.0
) y0 [ i1 ] = u0 [ i1 ] ; else y0 [ i1 ] =
Date0501_QuadModel_Thiago_acc_rt_TDelayInterpolate ( tMinusDelay , 0.0 , *
tBuffer , * uBuffer , iw_CircularBufSize [ i1 ] , & iw_Last [ i1 ] , iw_Tail
[ i1 ] , iw_Head [ i1 ] , _rtP -> P_26 , 0 , ( boolean_T ) (
ssIsMinorTimeStep ( S ) && ( ssGetTimeOfLastOutput ( S ) == ssGetT ( S ) ) )
) ; tBuffer ++ ; uBuffer ++ ; } } } ip43dgoh0e = _rtP -> P_27 * _rtB ->
az2wtq5aqx [ 0 ] ; cb55qal4hd = _rtP -> P_28 * _rtB -> az2wtq5aqx [ 1 ] ;
erv0ccpdoz = _rtP -> P_29 * _rtB -> az2wtq5aqx [ 2 ] ; h4fomgvj5n = _rtP ->
P_30 * _rtB -> az2wtq5aqx [ 3 ] ; a4qrcc3br4_idx_0 = ip43dgoh0e ;
a4qrcc3br4_idx_1 = cb55qal4hd ; a4qrcc3br4_idx_2 = erv0ccpdoz ;
a4qrcc3br4_idx_3 = h4fomgvj5n ; for ( i = 0 ; i < 4 ; i ++ ) { npg5zigbkh =
etztp0mvgy [ i + 12 ] * h4fomgvj5n + ( etztp0mvgy [ i + 8 ] * erv0ccpdoz + (
etztp0mvgy [ i + 4 ] * cb55qal4hd + etztp0mvgy [ i ] * ip43dgoh0e ) ) ;
hzpktek0ad [ i ] = npg5zigbkh ; } for ( i = 0 ; i < 4 ; i ++ ) { _rtB ->
eikp0wsyh4 [ i ] = 0.0 ; _rtB -> eikp0wsyh4 [ i ] += _rtB -> i51wezcj32 [ i ]
* hzpktek0ad [ 0 ] ; _rtB -> eikp0wsyh4 [ i ] += _rtB -> i51wezcj32 [ i + 4 ]
* hzpktek0ad [ 1 ] ; _rtB -> eikp0wsyh4 [ i ] += _rtB -> i51wezcj32 [ i + 8 ]
* hzpktek0ad [ 2 ] ; _rtB -> eikp0wsyh4 [ i ] += _rtB -> i51wezcj32 [ i + 12
] * hzpktek0ad [ 3 ] ; } } if ( ssIsSampleHit ( S , 1 , tid ) ) {
ssCallAccelRunBlock ( S , 2 , 41 , SS_CALL_MDL_OUTPUTS ) ;
ssCallAccelRunBlock ( S , 2 , 42 , SS_CALL_MDL_OUTPUTS ) ; } if (
ssIsContinuousTask ( S , tid ) ) { _rtB -> c0d553mc21 = muDoubleScalarCos (
gshcvxpyxd ) / 4.0 ; ip43dgoh0e = _rtB -> c0d553mc21 - _rtB -> crqrmutysw [ 0
] ; if ( ( _rtDW -> boh20pztcq >= ssGetT ( S ) ) && ( _rtDW -> lujb4z3nxo >=
ssGetT ( S ) ) ) { cb55qal4hd = 0.0 ; } else { nuudevirgz = _rtDW ->
boh20pztcq ; lastU = & _rtDW -> d2gontfkb3 ; if ( _rtDW -> boh20pztcq < _rtDW
-> lujb4z3nxo ) { if ( _rtDW -> lujb4z3nxo < ssGetT ( S ) ) { nuudevirgz =
_rtDW -> lujb4z3nxo ; lastU = & _rtDW -> b1luw3xjw0 ; } } else { if ( _rtDW
-> boh20pztcq >= ssGetT ( S ) ) { nuudevirgz = _rtDW -> lujb4z3nxo ; lastU =
& _rtDW -> b1luw3xjw0 ; } } cb55qal4hd = ( _rtB -> c0d553mc21 - * lastU ) / (
ssGetT ( S ) - nuudevirgz ) ; } cb55qal4hd = ( _rtB -> np4c234dzr [ 0 ] -
cb55qal4hd ) - _rtP -> P_31 * ip43dgoh0e ; erv0ccpdoz = ( ( ip43dgoh0e - (
cb55qal4hd + ip43dgoh0e ) ) - 0.1 * cb55qal4hd ) * ( 1.65 / _rtB ->
c4oqtnupxw ) ; _rtB -> o3ysyitf3k = muDoubleScalarSin ( gshcvxpyxd ) / 4.0 ;
ip43dgoh0e = _rtB -> o3ysyitf3k - _rtB -> crqrmutysw [ 1 ] ; if ( ( _rtDW ->
mdccoohgln >= ssGetT ( S ) ) && ( _rtDW -> gld3li1xxa >= ssGetT ( S ) ) ) {
cb55qal4hd = 0.0 ; } else { nuudevirgz = _rtDW -> mdccoohgln ; lastU = &
_rtDW -> f3cfqf4mix ; if ( _rtDW -> mdccoohgln < _rtDW -> gld3li1xxa ) { if (
_rtDW -> gld3li1xxa < ssGetT ( S ) ) { nuudevirgz = _rtDW -> gld3li1xxa ;
lastU = & _rtDW -> hbqvvmbhg4 ; } } else { if ( _rtDW -> mdccoohgln >= ssGetT
( S ) ) { nuudevirgz = _rtDW -> gld3li1xxa ; lastU = & _rtDW -> hbqvvmbhg4 ;
} } cb55qal4hd = ( _rtB -> o3ysyitf3k - * lastU ) / ( ssGetT ( S ) -
nuudevirgz ) ; } cb55qal4hd = ( _rtB -> np4c234dzr [ 1 ] - cb55qal4hd ) -
_rtP -> P_32 * ip43dgoh0e ; h4fomgvj5n = ( ( ip43dgoh0e - ( cb55qal4hd +
ip43dgoh0e ) ) - 0.1 * cb55qal4hd ) * ( 1.65 / _rtB -> c4oqtnupxw ) ;
ip43dgoh0e = muDoubleScalarCos ( _rtB -> fl4bn1wmqu [ 2 ] ) ; cb55qal4hd =
muDoubleScalarSin ( _rtB -> fl4bn1wmqu [ 2 ] ) ; nuudevirgz = cb55qal4hd *
erv0ccpdoz - ip43dgoh0e * h4fomgvj5n ; if ( nuudevirgz > 1.0 ) { nuudevirgz =
1.0 ; } else { if ( nuudevirgz < - 1.0 ) { nuudevirgz = - 1.0 ; } }
gmesnt5dc1 = muDoubleScalarAsin ( nuudevirgz ) ; erv0ccpdoz = ip43dgoh0e *
erv0ccpdoz + cb55qal4hd * h4fomgvj5n ; nuudevirgz = erv0ccpdoz /
muDoubleScalarCos ( gmesnt5dc1 ) ; if ( nuudevirgz > 1.0 ) { nuudevirgz = 1.0
; } else { if ( nuudevirgz < - 1.0 ) { nuudevirgz = - 1.0 ; } } _rtB ->
bgh02yljjm = muDoubleScalarAsin ( nuudevirgz ) ; gshcvxpyxd =
muDoubleScalarSin ( gshcvxpyxd ) ; if ( gmesnt5dc1 > _rtP -> P_33 ) { _rtB ->
f1sa0cbc35 [ 0 ] = _rtP -> P_33 ; } else if ( gmesnt5dc1 < _rtP -> P_34 ) {
_rtB -> f1sa0cbc35 [ 0 ] = _rtP -> P_34 ; } else { _rtB -> f1sa0cbc35 [ 0 ] =
gmesnt5dc1 ; } if ( _rtB -> bgh02yljjm > _rtP -> P_33 ) { _rtB -> f1sa0cbc35
[ 1 ] = _rtP -> P_33 ; } else if ( _rtB -> bgh02yljjm < _rtP -> P_34 ) { _rtB
-> f1sa0cbc35 [ 1 ] = _rtP -> P_34 ; } else { _rtB -> f1sa0cbc35 [ 1 ] = _rtB
-> bgh02yljjm ; } if ( gshcvxpyxd > _rtP -> P_33 ) { _rtB -> f1sa0cbc35 [ 2 ]
= _rtP -> P_33 ; } else if ( gshcvxpyxd < _rtP -> P_34 ) { _rtB -> f1sa0cbc35
[ 2 ] = _rtP -> P_34 ; } else { _rtB -> f1sa0cbc35 [ 2 ] = gshcvxpyxd ; } }
if ( ssIsSampleHit ( S , 1 , tid ) ) { ssCallAccelRunBlock ( S , 2 , 65 ,
SS_CALL_MDL_OUTPUTS ) ; ssCallAccelRunBlock ( S , 2 , 66 ,
SS_CALL_MDL_OUTPUTS ) ; ssCallAccelRunBlock ( S , 2 , 67 ,
SS_CALL_MDL_OUTPUTS ) ; ssCallAccelRunBlock ( S , 2 , 68 ,
SS_CALL_MDL_OUTPUTS ) ; ssCallAccelRunBlock ( S , 2 , 69 ,
SS_CALL_MDL_OUTPUTS ) ; ssCallAccelRunBlock ( S , 2 , 70 ,
SS_CALL_MDL_OUTPUTS ) ; } if ( ssIsContinuousTask ( S , tid ) ) { fegfldezgg
[ 0 ] = _rtX -> anaapvbrhg [ 0 ] ; fegfldezgg [ 1 ] = _rtX -> anaapvbrhg [ 1
] ; fegfldezgg [ 2 ] = _rtX -> anaapvbrhg [ 2 ] ; _rtB -> elfzxcrqe0 [ 0 ] =
_rtX -> anaapvbrhg [ 0 ] - _rtB -> maivwxyydn [ 0 ] ; _rtB -> elfzxcrqe0 [ 1
] = _rtX -> anaapvbrhg [ 1 ] - _rtB -> maivwxyydn [ 1 ] ; _rtB -> elfzxcrqe0
[ 2 ] = _rtX -> anaapvbrhg [ 2 ] - _rtB -> maivwxyydn [ 2 ] ; } if (
ssIsSampleHit ( S , 1 , tid ) ) { ssCallAccelRunBlock ( S , 2 , 74 ,
SS_CALL_MDL_OUTPUTS ) ; } if ( ssIsContinuousTask ( S , tid ) ) { _rtB ->
ck4hyyzuzm = ( ( _rtP -> P_39 * _rtB -> f1sa0cbc35 [ 0 ] - _rtB -> fl4bn1wmqu
[ 0 ] ) * _rtP -> P_40 - _rtX -> mxidsv51if ) * _rtP -> P_42 ; _rtB ->
h3aapfna3c = ( ( _rtP -> P_36 * _rtB -> f1sa0cbc35 [ 0 ] - _rtB -> fl4bn1wmqu
[ 0 ] ) * _rtP -> P_37 + _rtX -> d3b435ieqt ) + _rtB -> ck4hyyzuzm ; _rtB ->
fxbyzsjzri = ( ( _rtP -> P_46 * _rtB -> f1sa0cbc35 [ 1 ] - _rtB -> fl4bn1wmqu
[ 1 ] ) * _rtP -> P_47 - _rtX -> abtaamq0dn ) * _rtP -> P_49 ; _rtB ->
a32fh1fotb = ( ( _rtP -> P_43 * _rtB -> f1sa0cbc35 [ 1 ] - _rtB -> fl4bn1wmqu
[ 1 ] ) * _rtP -> P_44 + _rtX -> pmidr2fzsd ) + _rtB -> fxbyzsjzri ; _rtB ->
cycpvzgfbm = ( ( _rtP -> P_53 * _rtB -> f1sa0cbc35 [ 2 ] - _rtB -> fl4bn1wmqu
[ 2 ] ) * _rtP -> P_54 - _rtX -> olsaxgcylj ) * _rtP -> P_56 ; _rtB ->
ggylqlmbqp = ( ( _rtP -> P_50 * _rtB -> f1sa0cbc35 [ 2 ] - _rtB -> fl4bn1wmqu
[ 2 ] ) * _rtP -> P_51 + _rtX -> ljx4hjxhlm ) + _rtB -> cycpvzgfbm ; } if (
ssIsSampleHit ( S , 1 , tid ) ) { ssCallAccelRunBlock ( S , 2 , 108 ,
SS_CALL_MDL_OUTPUTS ) ; ssCallAccelRunBlock ( S , 2 , 109 ,
SS_CALL_MDL_OUTPUTS ) ; } if ( ssIsSampleHit ( S , 2 , tid ) ) { memcpy ( &
diznsy5jit [ 0 ] , & _rtP -> P_57 [ 0 ] , 9U * sizeof ( real_T ) ) ; } if (
ssIsContinuousTask ( S , tid ) && ssIsSpecialSampleHit ( S , 2 , 0 , tid ) )
{ _rtB -> eancbzmw2v [ 0 ] = _rtB -> elfzxcrqe0 [ 0 ] ; _rtB -> eancbzmw2v [
1 ] = _rtB -> elfzxcrqe0 [ 1 ] ; _rtB -> eancbzmw2v [ 2 ] = _rtB ->
elfzxcrqe0 [ 2 ] ; } if ( ssIsSampleHit ( S , 2 , tid ) ) { for ( i = 0 ; i <
3 ; i ++ ) { fvvgxc4onr [ i ] = diznsy5jit [ i + 6 ] * _rtB -> eancbzmw2v [ 2
] + ( diznsy5jit [ i + 3 ] * _rtB -> eancbzmw2v [ 1 ] + diznsy5jit [ i ] *
_rtB -> eancbzmw2v [ 0 ] ) ; } _rtB -> dq04431le5 [ 0 ] = _rtP -> P_58 *
fvvgxc4onr [ 0 ] ; _rtB -> dq04431le5 [ 1 ] = _rtP -> P_58 * fvvgxc4onr [ 1 ]
; _rtB -> dq04431le5 [ 2 ] = _rtP -> P_58 * fvvgxc4onr [ 2 ] ; } if (
ssIsSampleHit ( S , 1 , tid ) ) { if ( ssIsSpecialSampleHit ( S , 2 , 1 , tid
) ) { _rtB -> cjkcuifh31 [ 0 ] = _rtDW -> k0e0pysv1e [ 0 ] ; _rtB ->
cjkcuifh31 [ 1 ] = _rtDW -> k0e0pysv1e [ 1 ] ; _rtB -> cjkcuifh31 [ 2 ] =
_rtDW -> k0e0pysv1e [ 2 ] ; } if ( _rtB -> cjkcuifh31 [ 0 ] > _rtP -> P_60 )
{ _rtB -> oat10ddnmj [ 0 ] = _rtP -> P_60 ; } else if ( _rtB -> cjkcuifh31 [
0 ] < _rtP -> P_61 ) { _rtB -> oat10ddnmj [ 0 ] = _rtP -> P_61 ; } else {
_rtB -> oat10ddnmj [ 0 ] = _rtB -> cjkcuifh31 [ 0 ] ; } if ( _rtB ->
cjkcuifh31 [ 1 ] > _rtP -> P_60 ) { _rtB -> oat10ddnmj [ 1 ] = _rtP -> P_60 ;
} else if ( _rtB -> cjkcuifh31 [ 1 ] < _rtP -> P_61 ) { _rtB -> oat10ddnmj [
1 ] = _rtP -> P_61 ; } else { _rtB -> oat10ddnmj [ 1 ] = _rtB -> cjkcuifh31 [
1 ] ; } if ( _rtB -> cjkcuifh31 [ 2 ] > _rtP -> P_60 ) { _rtB -> oat10ddnmj [
2 ] = _rtP -> P_60 ; } else if ( _rtB -> cjkcuifh31 [ 2 ] < _rtP -> P_61 ) {
_rtB -> oat10ddnmj [ 2 ] = _rtP -> P_61 ; } else { _rtB -> oat10ddnmj [ 2 ] =
_rtB -> cjkcuifh31 [ 2 ] ; } memcpy ( & _rtB -> pm4j3jf0ja [ 0 ] , & _rtP ->
P_62 [ 0 ] , 9U * sizeof ( real_T ) ) ; } if ( ssIsContinuousTask ( S , tid )
) { ljqyvq454t [ 0 ] = _rtB -> h3aapfna3c ; ljqyvq454t [ 1 ] = _rtB ->
a32fh1fotb ; ljqyvq454t [ 2 ] = _rtB -> ggylqlmbqp ; for ( i = 0 ; i < 3 ; i
++ ) { aejymufvap [ i ] = _rtB -> pm4j3jf0ja [ i + 6 ] * _rtB -> ggylqlmbqp +
( _rtB -> pm4j3jf0ja [ i + 3 ] * _rtB -> a32fh1fotb + _rtB -> pm4j3jf0ja [ i
] * _rtB -> h3aapfna3c ) ; } for ( i = 0 ; i < 3 ; i ++ ) { _rtB ->
gq23ns2imo [ i ] = ( _rtB -> oat10ddnmj [ i ] - ( ( _rtB -> pm4j3jf0ja [ i +
3 ] * _rtB -> a32fh1fotb + _rtB -> pm4j3jf0ja [ i ] * _rtB -> h3aapfna3c ) +
_rtB -> pm4j3jf0ja [ i + 6 ] * _rtB -> ggylqlmbqp ) ) + gg2u4if1k1 [ i ] ; }
_rtB -> bo10f5adsd = ( _rtB -> f1sa0cbc35 [ 0 ] - _rtB -> fl4bn1wmqu [ 0 ] )
* _rtP -> P_63 ; _rtB -> pfqrkpg5dx = ( _rtB -> f1sa0cbc35 [ 1 ] - _rtB ->
fl4bn1wmqu [ 1 ] ) * _rtP -> P_64 ; _rtB -> fygzx32vw5 = ( _rtB -> f1sa0cbc35
[ 2 ] - _rtB -> fl4bn1wmqu [ 2 ] ) * _rtP -> P_65 ; } if ( ssIsSampleHit ( S
, 1 , tid ) ) { for ( i = 0 ; i < 9 ; i ++ ) { _rtB -> joeqwc11fy [ i ] =
_rtP -> P_66 [ i ] ; _rtB -> gr5my0drs2 [ i ] = _rtP -> P_67 [ i ] ; } } if (
ssIsContinuousTask ( S , tid ) ) { for ( i = 0 ; i < 3 ; i ++ ) { aejymufvap
[ i ] = _rtB -> joeqwc11fy [ i + 6 ] * fegfldezgg [ 2 ] + ( _rtB ->
joeqwc11fy [ i + 3 ] * fegfldezgg [ 1 ] + _rtB -> joeqwc11fy [ i ] *
fegfldezgg [ 0 ] ) ; } gg2u4if1k1 [ 0 ] += _rtB -> oat10ddnmj [ 0 ] ;
gg2u4if1k1 [ 1 ] += _rtB -> oat10ddnmj [ 1 ] ; npg5zigbkh = gg2u4if1k1 [ 2 ]
+ _rtB -> oat10ddnmj [ 2 ] ; for ( i = 0 ; i < 3 ; i ++ ) { ljqyvq454t [ i ]
= _rtB -> gr5my0drs2 [ i + 6 ] * npg5zigbkh + ( _rtB -> gr5my0drs2 [ i + 3 ]
* gg2u4if1k1 [ 1 ] + _rtB -> gr5my0drs2 [ i ] * gg2u4if1k1 [ 0 ] ) ; } for (
i = 0 ; i < 3 ; i ++ ) { tmp [ i ] = _rtB -> gr5my0drs2 [ i + 6 ] *
npg5zigbkh + ( _rtB -> gr5my0drs2 [ i + 3 ] * gg2u4if1k1 [ 1 ] + _rtB ->
gr5my0drs2 [ i ] * gg2u4if1k1 [ 0 ] ) ; } for ( i = 0 ; i < 3 ; i ++ ) {
fvvgxc4onr [ i ] = _rtB -> joeqwc11fy [ i + 6 ] * fegfldezgg [ 2 ] + ( _rtB
-> joeqwc11fy [ i + 3 ] * fegfldezgg [ 1 ] + _rtB -> joeqwc11fy [ i ] *
fegfldezgg [ 0 ] ) ; } _rtB -> gu4qmbqny2 [ 0 ] = tmp [ 0 ] + fvvgxc4onr [ 0
] ; _rtB -> gu4qmbqny2 [ 1 ] = tmp [ 1 ] + fvvgxc4onr [ 1 ] ; _rtB ->
gu4qmbqny2 [ 2 ] = tmp [ 2 ] + fvvgxc4onr [ 2 ] ; nuudevirgz = erv0ccpdoz /
muDoubleScalarCos ( gmesnt5dc1 ) ; if ( nuudevirgz > 1.0 ) { nuudevirgz = 1.0
; } else { if ( nuudevirgz < - 1.0 ) { nuudevirgz = - 1.0 ; } } _rtB ->
aa2n4x0fzs = muDoubleScalarAsin ( nuudevirgz ) ; } if ( ssIsSampleHit ( S , 1
, tid ) ) { ssCallAccelRunBlock ( S , 2 , 133 , SS_CALL_MDL_OUTPUTS ) ; } if
( ssIsContinuousTask ( S , tid ) ) { cb55qal4hd = _rtP -> P_73 * _rtX ->
lclvofzdkw ; if ( ssGetTaskTime ( S , 0 ) < _rtP -> P_74 ) { _rtB ->
i4h5t5xorw = _rtP -> P_75 ; } else { _rtB -> i4h5t5xorw = _rtP -> P_76 ; } if
( ( _rtDW -> jgbpzk0jm2 >= ssGetT ( S ) ) && ( _rtDW -> aesgfko0p3 >= ssGetT
( S ) ) ) { gshcvxpyxd = 0.0 ; } else { nuudevirgz = _rtDW -> jgbpzk0jm2 ;
lastU = & _rtDW -> djexgjadxg ; if ( _rtDW -> jgbpzk0jm2 < _rtDW ->
aesgfko0p3 ) { if ( _rtDW -> aesgfko0p3 < ssGetT ( S ) ) { nuudevirgz = _rtDW
-> aesgfko0p3 ; lastU = & _rtDW -> oxiu4gui43 ; } } else { if ( _rtDW ->
jgbpzk0jm2 >= ssGetT ( S ) ) { nuudevirgz = _rtDW -> aesgfko0p3 ; lastU = &
_rtDW -> oxiu4gui43 ; } } gshcvxpyxd = ( _rtB -> i4h5t5xorw - * lastU ) / (
ssGetT ( S ) - nuudevirgz ) ; } if ( gshcvxpyxd > _rtP -> P_77 ) { gshcvxpyxd
= _rtP -> P_77 ; } else { if ( gshcvxpyxd < _rtP -> P_78 ) { gshcvxpyxd =
_rtP -> P_78 ; } } a4qrcc3br4_idx_0 = ( _rtP -> P_69 * _rtX -> ocryheeszg +
_rtB -> eikp0wsyh4 [ 0 ] ) + gshcvxpyxd ; a4qrcc3br4_idx_1 = ( _rtP -> P_71 *
_rtX -> byywf5sz3i + _rtB -> eikp0wsyh4 [ 1 ] ) + 0.0 ; a4qrcc3br4_idx_2 = (
_rtB -> eikp0wsyh4 [ 2 ] + cb55qal4hd ) + gshcvxpyxd ; a4qrcc3br4_idx_3 = (
_rtB -> eikp0wsyh4 [ 3 ] + cb55qal4hd ) + 0.0 ; } if ( ssIsSampleHit ( S , 1
, tid ) ) { adoq34mqbb = _rtP -> P_82 * _rtDW -> g05avm0rip ; iz2k5ya5ta =
_rtP -> P_86 * _rtDW -> cclupfu4ez ; inyi50y04a = _rtP -> P_90 * _rtDW ->
ojclvvbjqu ; memcpy ( & _rtB -> iu2532friw [ 0 ] , & _rtP -> P_91 [ 0 ] , 9U
* sizeof ( real_T ) ) ; } if ( ssIsContinuousTask ( S , tid ) ) { for ( i = 0
; i < 3 ; i ++ ) { aejymufvap [ i ] = _rtB -> iu2532friw [ i + 6 ] * _rtB ->
maivwxyydn [ 2 ] + ( _rtB -> iu2532friw [ i + 3 ] * _rtB -> maivwxyydn [ 1 ]
+ _rtB -> iu2532friw [ i ] * _rtB -> maivwxyydn [ 0 ] ) ; } for ( i = 0 ; i <
3 ; i ++ ) { tmp [ i ] = _rtB -> iu2532friw [ i + 6 ] * _rtB -> maivwxyydn [
2 ] + ( _rtB -> iu2532friw [ i + 3 ] * _rtB -> maivwxyydn [ 1 ] + _rtB ->
iu2532friw [ i ] * _rtB -> maivwxyydn [ 0 ] ) ; } gshcvxpyxd = _rtB ->
maivwxyydn [ 0 ] * tmp [ 2 ] ; npg5zigbkh = aejymufvap [ 0 ] ; cb55qal4hd =
aejymufvap [ 1 ] ; nuudevirgz = aejymufvap [ 0 ] ; aejymufvap [ 0 ] = _rtB ->
maivwxyydn [ 1 ] * aejymufvap [ 2 ] - _rtB -> maivwxyydn [ 2 ] * aejymufvap [
1 ] ; aejymufvap [ 0 ] = _rtP -> P_92 * a4qrcc3br4_idx_1 - aejymufvap [ 0 ] ;
aejymufvap [ 1 ] = _rtP -> P_93 * a4qrcc3br4_idx_2 - ( _rtB -> maivwxyydn [ 2
] * npg5zigbkh - gshcvxpyxd ) ; aejymufvap [ 2 ] = a4qrcc3br4_idx_3 - ( _rtB
-> maivwxyydn [ 0 ] * cb55qal4hd - _rtB -> maivwxyydn [ 1 ] * nuudevirgz ) ;
} if ( ssIsSampleHit ( S , 1 , tid ) ) { memcpy ( & _rtB -> g3fwejnijf [ 0 ]
, & _rtP -> P_94 [ 0 ] , 9U * sizeof ( real_T ) ) ; } if ( ssIsContinuousTask
( S , tid ) ) { for ( i = 0 ; i < 3 ; i ++ ) { ljqyvq454t [ i ] = _rtB ->
g3fwejnijf [ i + 6 ] * aejymufvap [ 2 ] + ( _rtB -> g3fwejnijf [ i + 3 ] *
aejymufvap [ 1 ] + _rtB -> g3fwejnijf [ i ] * aejymufvap [ 0 ] ) ; } } if (
ssIsSampleHit ( S , 1 , tid ) ) { _rtB -> o1znxitunf [ 0 ] = _rtP -> P_95 [ 0
] ; _rtB -> o1znxitunf [ 1 ] = _rtP -> P_95 [ 1 ] ; _rtB -> o1znxitunf [ 2 ]
= _rtP -> P_95 [ 2 ] ; _rtB -> o1znxitunf [ 3 ] = _rtP -> P_95 [ 3 ] ; } if (
ssIsContinuousTask ( S , tid ) ) { hzpktek0ad [ 0 ] *= _rtB -> o1znxitunf [ 0
] ; hzpktek0ad [ 1 ] *= _rtB -> o1znxitunf [ 1 ] ; hzpktek0ad [ 2 ] *= _rtB
-> o1znxitunf [ 2 ] ; gmesnt5dc1 = ( ( hzpktek0ad [ 0 ] + hzpktek0ad [ 1 ] )
+ hzpktek0ad [ 2 ] ) + _rtB -> o1znxitunf [ 3 ] * hzpktek0ad [ 3 ] ; _rtB ->
kngt2l52e3 [ 0 ] = _rtP -> P_96 * gmesnt5dc1 * _rtB -> maivwxyydn [ 1 ] +
ljqyvq454t [ 0 ] ; _rtB -> kngt2l52e3 [ 1 ] = _rtP -> P_97 * gmesnt5dc1 *
_rtB -> maivwxyydn [ 0 ] + ljqyvq454t [ 1 ] ; _rtB -> kngt2l52e3 [ 2 ] =
ljqyvq454t [ 2 ] + 0.0 ; } if ( ssIsSampleHit ( S , 1 , tid ) ) { _rtB ->
mm53i31ky2 = _rtP -> P_98 * iz2k5ya5ta ; _rtB -> mrol05thcy = _rtP -> P_99 *
adoq34mqbb ; _rtB -> kps5fe005v = _rtP -> P_100 * inyi50y04a ; } if (
ssIsContinuousTask ( S , tid ) ) { gmesnt5dc1 = muDoubleScalarCos ( _rtB ->
fl4bn1wmqu [ 0 ] ) ; ip43dgoh0e = muDoubleScalarSin ( _rtB -> fl4bn1wmqu [ 1
] ) ; cb55qal4hd = muDoubleScalarCos ( _rtB -> fl4bn1wmqu [ 2 ] ) ;
gshcvxpyxd = muDoubleScalarSin ( _rtB -> fl4bn1wmqu [ 0 ] ) ; erv0ccpdoz =
muDoubleScalarSin ( _rtB -> fl4bn1wmqu [ 2 ] ) ; _rtB -> efrpdb2llp = (
gmesnt5dc1 * ip43dgoh0e * cb55qal4hd + gshcvxpyxd * erv0ccpdoz ) *
a4qrcc3br4_idx_0 * _rtP -> P_101 ; _rtB -> he10hnsld4 = ( gmesnt5dc1 *
ip43dgoh0e * erv0ccpdoz - gshcvxpyxd * cb55qal4hd ) * a4qrcc3br4_idx_0 * _rtP
-> P_102 ; nh5etxksnt = gmesnt5dc1 * muDoubleScalarCos ( _rtB -> fl4bn1wmqu [
1 ] ) * a4qrcc3br4_idx_0 * _rtP -> P_103 ; } if ( ssIsSampleHit ( S , 1 , tid
) ) { _rtB -> nvxvftase4 = _rtP -> P_104 ; } if ( ssIsContinuousTask ( S ,
tid ) ) { _rtB -> pzzpvzwlyt = _rtB -> nvxvftase4 + nh5etxksnt ; } if (
ssIsSampleHit ( S , 1 , tid ) ) { ssCallAccelRunBlock ( S , 2 , 193 ,
SS_CALL_MDL_OUTPUTS ) ; ssCallAccelRunBlock ( S , 2 , 194 ,
SS_CALL_MDL_OUTPUTS ) ; } UNUSED_PARAMETER ( tid ) ; }
#define MDL_UPDATE
static void mdlUpdate ( SimStruct * S , int_T tid ) { real_T * lastU ; char_T
* sErr ; gsegsuzktn * _rtB ; pcuunorfmd * _rtP ; pfm4mif5gq * _rtDW ; _rtDW =
( ( pfm4mif5gq * ) ssGetRootDWork ( S ) ) ; _rtP = ( ( pcuunorfmd * )
ssGetDefaultParam ( S ) ) ; _rtB = ( ( gsegsuzktn * ) _ssGetBlockIO ( S ) ) ;
if ( ssIsContinuousTask ( S , tid ) ) { if ( _rtDW -> kmpk31urxn == ( rtInf )
) { _rtDW -> kmpk31urxn = ssGetT ( S ) ; lastU = & _rtDW -> hphjcnocig ; }
else if ( _rtDW -> ogdw1xlqmp == ( rtInf ) ) { _rtDW -> ogdw1xlqmp = ssGetT (
S ) ; lastU = & _rtDW -> jjecbj5saz ; } else if ( _rtDW -> kmpk31urxn < _rtDW
-> ogdw1xlqmp ) { _rtDW -> kmpk31urxn = ssGetT ( S ) ; lastU = & _rtDW ->
hphjcnocig ; } else { _rtDW -> ogdw1xlqmp = ssGetT ( S ) ; lastU = & _rtDW ->
jjecbj5saz ; } * lastU = _rtB -> b3px43da1a ; { real_T * * uBuffer = ( real_T
* * ) & _rtDW -> nulonpacio . TUbufferPtrs [ 0 ] ; real_T * * tBuffer = (
real_T * * ) & _rtDW -> nulonpacio . TUbufferPtrs [ 4 ] ; real_T simTime =
ssGetT ( S ) ; _rtDW -> ijfvsgcbuh . Head [ 0 ] = ( ( _rtDW -> ijfvsgcbuh .
Head [ 0 ] < ( _rtDW -> ijfvsgcbuh . CircularBufSize [ 0 ] - 1 ) ) ? ( _rtDW
-> ijfvsgcbuh . Head [ 0 ] + 1 ) : 0 ) ; if ( _rtDW -> ijfvsgcbuh . Head [ 0
] == _rtDW -> ijfvsgcbuh . Tail [ 0 ] ) { if ( !
Date0501_QuadModel_Thiago_acc_rt_TDelayUpdateTailOrGrowBuf ( & _rtDW ->
ijfvsgcbuh . CircularBufSize [ 0 ] , & _rtDW -> ijfvsgcbuh . Tail [ 0 ] , &
_rtDW -> ijfvsgcbuh . Head [ 0 ] , & _rtDW -> ijfvsgcbuh . Last [ 0 ] ,
simTime - _rtP -> P_25 , tBuffer , uBuffer , ( NULL ) , ( boolean_T ) 0 ,
false , & _rtDW -> ijfvsgcbuh . MaxNewBufSize ) ) { ssSetErrorStatus ( S ,
"tdelay memory allocation error" ) ; return ; } } ( * tBuffer ++ ) [ _rtDW ->
ijfvsgcbuh . Head [ 0 ] ] = simTime ; ( * uBuffer ++ ) [ _rtDW -> ijfvsgcbuh
. Head [ 0 ] ] = _rtB -> js2p3cmsfo [ 0 ] ; _rtDW -> ijfvsgcbuh . Head [ 1 ]
= ( ( _rtDW -> ijfvsgcbuh . Head [ 1 ] < ( _rtDW -> ijfvsgcbuh .
CircularBufSize [ 1 ] - 1 ) ) ? ( _rtDW -> ijfvsgcbuh . Head [ 1 ] + 1 ) : 0
) ; if ( _rtDW -> ijfvsgcbuh . Head [ 1 ] == _rtDW -> ijfvsgcbuh . Tail [ 1 ]
) { if ( ! Date0501_QuadModel_Thiago_acc_rt_TDelayUpdateTailOrGrowBuf ( &
_rtDW -> ijfvsgcbuh . CircularBufSize [ 1 ] , & _rtDW -> ijfvsgcbuh . Tail [
1 ] , & _rtDW -> ijfvsgcbuh . Head [ 1 ] , & _rtDW -> ijfvsgcbuh . Last [ 1 ]
, simTime - _rtP -> P_25 , tBuffer , uBuffer , ( NULL ) , ( boolean_T ) 0 ,
false , & _rtDW -> ijfvsgcbuh . MaxNewBufSize ) ) { ssSetErrorStatus ( S ,
"tdelay memory allocation error" ) ; return ; } } ( * tBuffer ++ ) [ _rtDW ->
ijfvsgcbuh . Head [ 1 ] ] = simTime ; ( * uBuffer ++ ) [ _rtDW -> ijfvsgcbuh
. Head [ 1 ] ] = _rtB -> js2p3cmsfo [ 1 ] ; _rtDW -> ijfvsgcbuh . Head [ 2 ]
= ( ( _rtDW -> ijfvsgcbuh . Head [ 2 ] < ( _rtDW -> ijfvsgcbuh .
CircularBufSize [ 2 ] - 1 ) ) ? ( _rtDW -> ijfvsgcbuh . Head [ 2 ] + 1 ) : 0
) ; if ( _rtDW -> ijfvsgcbuh . Head [ 2 ] == _rtDW -> ijfvsgcbuh . Tail [ 2 ]
) { if ( ! Date0501_QuadModel_Thiago_acc_rt_TDelayUpdateTailOrGrowBuf ( &
_rtDW -> ijfvsgcbuh . CircularBufSize [ 2 ] , & _rtDW -> ijfvsgcbuh . Tail [
2 ] , & _rtDW -> ijfvsgcbuh . Head [ 2 ] , & _rtDW -> ijfvsgcbuh . Last [ 2 ]
, simTime - _rtP -> P_25 , tBuffer , uBuffer , ( NULL ) , ( boolean_T ) 0 ,
false , & _rtDW -> ijfvsgcbuh . MaxNewBufSize ) ) { ssSetErrorStatus ( S ,
"tdelay memory allocation error" ) ; return ; } } ( * tBuffer ++ ) [ _rtDW ->
ijfvsgcbuh . Head [ 2 ] ] = simTime ; ( * uBuffer ++ ) [ _rtDW -> ijfvsgcbuh
. Head [ 2 ] ] = _rtB -> js2p3cmsfo [ 2 ] ; _rtDW -> ijfvsgcbuh . Head [ 3 ]
= ( ( _rtDW -> ijfvsgcbuh . Head [ 3 ] < ( _rtDW -> ijfvsgcbuh .
CircularBufSize [ 3 ] - 1 ) ) ? ( _rtDW -> ijfvsgcbuh . Head [ 3 ] + 1 ) : 0
) ; if ( _rtDW -> ijfvsgcbuh . Head [ 3 ] == _rtDW -> ijfvsgcbuh . Tail [ 3 ]
) { if ( ! Date0501_QuadModel_Thiago_acc_rt_TDelayUpdateTailOrGrowBuf ( &
_rtDW -> ijfvsgcbuh . CircularBufSize [ 3 ] , & _rtDW -> ijfvsgcbuh . Tail [
3 ] , & _rtDW -> ijfvsgcbuh . Head [ 3 ] , & _rtDW -> ijfvsgcbuh . Last [ 3 ]
, simTime - _rtP -> P_25 , tBuffer , uBuffer , ( NULL ) , ( boolean_T ) 0 ,
false , & _rtDW -> ijfvsgcbuh . MaxNewBufSize ) ) { ssSetErrorStatus ( S ,
"tdelay memory allocation error" ) ; return ; } } ( * tBuffer ) [ _rtDW ->
ijfvsgcbuh . Head [ 3 ] ] = simTime ; ( * uBuffer ) [ _rtDW -> ijfvsgcbuh .
Head [ 3 ] ] = _rtB -> js2p3cmsfo [ 3 ] ; } } if ( ssIsContinuousTask ( S ,
tid ) ) { if ( _rtDW -> boh20pztcq == ( rtInf ) ) { _rtDW -> boh20pztcq =
ssGetT ( S ) ; lastU = & _rtDW -> d2gontfkb3 ; } else if ( _rtDW ->
lujb4z3nxo == ( rtInf ) ) { _rtDW -> lujb4z3nxo = ssGetT ( S ) ; lastU = &
_rtDW -> b1luw3xjw0 ; } else if ( _rtDW -> boh20pztcq < _rtDW -> lujb4z3nxo )
{ _rtDW -> boh20pztcq = ssGetT ( S ) ; lastU = & _rtDW -> d2gontfkb3 ; } else
{ _rtDW -> lujb4z3nxo = ssGetT ( S ) ; lastU = & _rtDW -> b1luw3xjw0 ; } *
lastU = _rtB -> c0d553mc21 ; if ( _rtDW -> mdccoohgln == ( rtInf ) ) { _rtDW
-> mdccoohgln = ssGetT ( S ) ; lastU = & _rtDW -> f3cfqf4mix ; } else if (
_rtDW -> gld3li1xxa == ( rtInf ) ) { _rtDW -> gld3li1xxa = ssGetT ( S ) ;
lastU = & _rtDW -> hbqvvmbhg4 ; } else if ( _rtDW -> mdccoohgln < _rtDW ->
gld3li1xxa ) { _rtDW -> mdccoohgln = ssGetT ( S ) ; lastU = & _rtDW ->
f3cfqf4mix ; } else { _rtDW -> gld3li1xxa = ssGetT ( S ) ; lastU = & _rtDW ->
hbqvvmbhg4 ; } * lastU = _rtB -> o3ysyitf3k ; } if ( ssIsSampleHit ( S , 1 ,
tid ) ) { sErr = GetErrorBuffer ( & _rtDW -> bgqcrmdnhx [ 0U ] ) ;
LibUpdate_Network ( & _rtDW -> bgqcrmdnhx [ 0U ] , & _rtB -> opfwoldhlb [ 0U
] , 48 ) ; if ( * sErr != 0 ) { ssSetErrorStatus ( S , sErr ) ;
ssSetStopRequested ( S , 1 ) ; } } if ( ssIsSampleHit ( S , 2 , tid ) ) {
_rtDW -> k0e0pysv1e [ 0 ] = _rtB -> dq04431le5 [ 0 ] ; _rtDW -> k0e0pysv1e [
1 ] = _rtB -> dq04431le5 [ 1 ] ; _rtDW -> k0e0pysv1e [ 2 ] = _rtB ->
dq04431le5 [ 2 ] ; } if ( ssIsContinuousTask ( S , tid ) ) { if ( _rtDW ->
jgbpzk0jm2 == ( rtInf ) ) { _rtDW -> jgbpzk0jm2 = ssGetT ( S ) ; lastU = &
_rtDW -> djexgjadxg ; } else if ( _rtDW -> aesgfko0p3 == ( rtInf ) ) { _rtDW
-> aesgfko0p3 = ssGetT ( S ) ; lastU = & _rtDW -> oxiu4gui43 ; } else if (
_rtDW -> jgbpzk0jm2 < _rtDW -> aesgfko0p3 ) { _rtDW -> jgbpzk0jm2 = ssGetT (
S ) ; lastU = & _rtDW -> djexgjadxg ; } else { _rtDW -> aesgfko0p3 = ssGetT (
S ) ; lastU = & _rtDW -> oxiu4gui43 ; } * lastU = _rtB -> i4h5t5xorw ; } if (
ssIsSampleHit ( S , 1 , tid ) ) { _rtDW -> g05avm0rip =
rt_nrand_Upu32_Yd_f_pw_snf ( & _rtDW -> ou3qtdg5gq ) * _rtP -> P_80 + _rtP ->
P_79 ; _rtDW -> cclupfu4ez = rt_nrand_Upu32_Yd_f_pw_snf ( & _rtDW ->
p4t1eiv534 ) * _rtP -> P_84 + _rtP -> P_83 ; _rtDW -> ojclvvbjqu =
rt_nrand_Upu32_Yd_f_pw_snf ( & _rtDW -> hsgdfd0evk ) * _rtP -> P_88 + _rtP ->
P_87 ; } UNUSED_PARAMETER ( tid ) ; }
#define MDL_DERIVATIVES
static void mdlDerivatives ( SimStruct * S ) { gsegsuzktn * _rtB ; pcuunorfmd
* _rtP ; pvlbxo1noz * _rtX ; dagv4455as * _rtXdot ; _rtXdot = ( ( dagv4455as
* ) ssGetdX ( S ) ) ; _rtX = ( ( pvlbxo1noz * ) ssGetContStates ( S ) ) ;
_rtP = ( ( pcuunorfmd * ) ssGetDefaultParam ( S ) ) ; _rtB = ( ( gsegsuzktn *
) _ssGetBlockIO ( S ) ) ; _rtXdot -> mdqhadbmzi [ 0 ] = _rtB -> kngt2l52e3 [
0 ] ; _rtXdot -> mdqhadbmzi [ 1 ] = _rtB -> kngt2l52e3 [ 1 ] ; _rtXdot ->
mdqhadbmzi [ 2 ] = _rtB -> kngt2l52e3 [ 2 ] ; _rtXdot -> kswdq12pcl = 0.0 ;
_rtXdot -> kswdq12pcl += _rtP -> P_4 * _rtX -> kswdq12pcl ; _rtXdot ->
kswdq12pcl += _rtB -> gq23ns2imo [ 0 ] ; _rtXdot -> bpy1pus01d = 0.0 ;
_rtXdot -> bpy1pus01d += _rtP -> P_6 * _rtX -> bpy1pus01d ; _rtXdot ->
bpy1pus01d += _rtB -> gq23ns2imo [ 1 ] ; _rtXdot -> pqy1xvwazu = 0.0 ;
_rtXdot -> pqy1xvwazu += _rtP -> P_8 * _rtX -> pqy1xvwazu ; _rtXdot ->
pqy1xvwazu += _rtB -> gq23ns2imo [ 2 ] ; _rtXdot -> dwor1dwwuw [ 0 ] = _rtB
-> maivwxyydn [ 0 ] ; _rtXdot -> dwor1dwwuw [ 1 ] = _rtB -> maivwxyydn [ 1 ]
; _rtXdot -> dwor1dwwuw [ 2 ] = _rtB -> maivwxyydn [ 2 ] ; _rtXdot ->
mseabyslip [ 0 ] = _rtB -> np4c234dzr [ 0 ] ; _rtXdot -> mseabyslip [ 1 ] =
_rtB -> np4c234dzr [ 1 ] ; _rtXdot -> mseabyslip [ 2 ] = _rtB -> np4c234dzr [
2 ] ; _rtXdot -> mfpumpgq5t [ 0 ] = _rtB -> efrpdb2llp ; _rtXdot ->
mfpumpgq5t [ 1 ] = _rtB -> he10hnsld4 ; _rtXdot -> mfpumpgq5t [ 2 ] = _rtB ->
pzzpvzwlyt ; _rtXdot -> anaapvbrhg [ 0 ] = _rtB -> gu4qmbqny2 [ 0 ] ; _rtXdot
-> anaapvbrhg [ 1 ] = _rtB -> gu4qmbqny2 [ 1 ] ; _rtXdot -> anaapvbrhg [ 2 ]
= _rtB -> gu4qmbqny2 [ 2 ] ; _rtXdot -> d3b435ieqt = _rtB -> bo10f5adsd ;
_rtXdot -> mxidsv51if = _rtB -> ck4hyyzuzm ; _rtXdot -> pmidr2fzsd = _rtB ->
pfqrkpg5dx ; _rtXdot -> abtaamq0dn = _rtB -> fxbyzsjzri ; _rtXdot ->
ljx4hjxhlm = _rtB -> fygzx32vw5 ; _rtXdot -> olsaxgcylj = _rtB -> cycpvzgfbm
; _rtXdot -> ocryheeszg = 0.0 ; _rtXdot -> ocryheeszg += _rtP -> P_68 * _rtX
-> ocryheeszg ; _rtXdot -> ocryheeszg += _rtB -> mrol05thcy ; _rtXdot ->
byywf5sz3i = 0.0 ; _rtXdot -> byywf5sz3i += _rtP -> P_70 * _rtX -> byywf5sz3i
; _rtXdot -> byywf5sz3i += _rtB -> mm53i31ky2 ; _rtXdot -> lclvofzdkw = 0.0 ;
_rtXdot -> lclvofzdkw += _rtP -> P_72 * _rtX -> lclvofzdkw ; _rtXdot ->
lclvofzdkw += _rtB -> kps5fe005v ; } static void mdlInitializeSizes (
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
pfm4mif5gq ) ) { ssSetErrorStatus ( S ,
"Unexpected error: Internal DWork sizes do "
"not match for accelerator mex file." ) ; } if ( ssGetSizeofGlobalBlockIO ( S
) != sizeof ( gsegsuzktn ) ) { ssSetErrorStatus ( S ,
"Unexpected error: Internal BlockIO sizes do "
"not match for accelerator mex file." ) ; } { int ssSizeofParams ;
ssGetSizeofParams ( S , & ssSizeofParams ) ; if ( ssSizeofParams != sizeof (
pcuunorfmd ) ) { static char msg [ 256 ] ; sprintf ( msg ,
"Unexpected error: Internal Parameters sizes do "
"not match for accelerator mex file." ) ; } } _ssSetDefaultParam ( S , (
real_T * ) & lgooc43qje ) ; rt_InitInfAndNaN ( sizeof ( real_T ) ) ; } static
void mdlInitializeSampleTimes ( SimStruct * S ) { } static void mdlTerminate
( SimStruct * S ) { }
#include "simulink.c"
