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
real_T emusy1jcfj ; real_T * lastU ; real_T lastTime ; real_T fegfldezgg [ 3
] ; real_T gg2u4if1k1 [ 3 ] ; real_T pvstwylnp5 ; real_T erv0ccpdoz ; real_T
ieaqk3zo1a ; real_T ip43dgoh0e ; real_T hzpktek0ad [ 4 ] ; real_T gmesnt5dc1
; real_T ljqyvq454t [ 3 ] ; real_T aejymufvap [ 3 ] ; real_T nh5etxksnt ;
real_T etztp0mvgy [ 16 ] ; real_T diznsy5jit [ 9 ] ; int32_T i ; real_T
o2vmpwsi0u [ 3 ] ; real_T tmp [ 3 ] ; real_T no2x3xsfnt_idx_0 ; real_T
no2x3xsfnt_idx_1 ; real_T no2x3xsfnt_idx_2 ; real_T no2x3xsfnt_idx_3 ; real_T
u0 ; pllwawth2v * _rtB ; lfrdjjasqz * _rtP ; isuplyng3d * _rtX ; fw4wgdftov *
_rtDW ; _rtDW = ( ( fw4wgdftov * ) ssGetRootDWork ( S ) ) ; _rtX = ( (
isuplyng3d * ) ssGetContStates ( S ) ) ; _rtP = ( ( lfrdjjasqz * )
ssGetDefaultParam ( S ) ) ; _rtB = ( ( pllwawth2v * ) _ssGetBlockIO ( S ) ) ;
if ( ssIsSampleHit ( S , 1 , tid ) ) { memcpy ( & _rtB -> nr4coujhhw [ 0 ] ,
& _rtP -> P_0 [ 0 ] , 9U * sizeof ( real_T ) ) ; } if ( ssIsContinuousTask (
S , tid ) ) { _rtB -> mc53goiddk [ 0 ] = _rtX -> iidw152zl0 [ 0 ] ; _rtB ->
mc53goiddk [ 1 ] = _rtX -> iidw152zl0 [ 1 ] ; _rtB -> mc53goiddk [ 2 ] = _rtX
-> iidw152zl0 [ 2 ] ; for ( i = 0 ; i < 3 ; i ++ ) { fegfldezgg [ i ] = _rtB
-> nr4coujhhw [ i + 6 ] * _rtB -> mc53goiddk [ 2 ] + ( _rtB -> nr4coujhhw [ i
+ 3 ] * _rtB -> mc53goiddk [ 1 ] + _rtB -> nr4coujhhw [ i ] * _rtB ->
mc53goiddk [ 0 ] ) ; } emusy1jcfj = 0.0 ; emusy1jcfj += _rtP -> P_3 * _rtX ->
awllp4kjam ; erv0ccpdoz = _rtP -> P_5 * _rtX -> g11z1syag2 ; pvstwylnp5 =
_rtP -> P_7 * _rtX -> bwqacranye ; gg2u4if1k1 [ 0 ] = _rtP -> P_8 *
emusy1jcfj * _rtP -> P_9 ; gg2u4if1k1 [ 1 ] = _rtP -> P_8 * erv0ccpdoz * _rtP
-> P_9 ; gg2u4if1k1 [ 2 ] = _rtP -> P_8 * pvstwylnp5 * _rtP -> P_9 ; for ( i
= 0 ; i < 3 ; i ++ ) { _rtB -> dy0jts3w2k [ i ] = ( ( _rtB -> nr4coujhhw [ i
+ 3 ] * _rtB -> mc53goiddk [ 1 ] + _rtB -> nr4coujhhw [ i ] * _rtB ->
mc53goiddk [ 0 ] ) + _rtB -> nr4coujhhw [ i + 6 ] * _rtB -> mc53goiddk [ 2 ]
) + gg2u4if1k1 [ i ] ; } } if ( ssIsSampleHit ( S , 1 , tid ) ) { for ( i = 0
; i < 16 ; i ++ ) { _rtB -> n3bbgc0hhm [ i ] = _rtP -> P_10 [ i ] ; _rtB ->
fk3ehp0ppr [ i ] = _rtP -> P_11 [ i ] ; } } if ( ssIsContinuousTask ( S , tid
) ) { if ( ssGetTaskTime ( S , 0 ) < _rtP -> P_12 ) { pvstwylnp5 = _rtP ->
P_13 ; } else { pvstwylnp5 = _rtP -> P_14 ; } } if ( ssIsSampleHit ( S , 1 ,
tid ) ) { memcpy ( & _rtB -> fwkmvd2sxh [ 0 ] , & _rtP -> P_15 [ 0 ] , sizeof
( real_T ) << 4U ) ; } if ( ssIsContinuousTask ( S , tid ) ) { for ( i = 0 ;
i < 16 ; i ++ ) { if ( pvstwylnp5 >= _rtP -> P_16 ) { etztp0mvgy [ i ] = _rtB
-> fk3ehp0ppr [ i ] ; } else { etztp0mvgy [ i ] = _rtB -> fwkmvd2sxh [ i ] ;
} } } if ( ssIsSampleHit ( S , 1 , tid ) ) { memcpy ( & _rtB -> agmy4ow3om [
0 ] , & _rtP -> P_17 [ 0 ] , sizeof ( real_T ) << 4U ) ; } if (
ssIsContinuousTask ( S , tid ) ) { _rtB -> o4ip5fjwmd [ 0 ] = _rtX ->
chvazs3t52 [ 0 ] ; _rtB -> o4ip5fjwmd [ 1 ] = _rtX -> chvazs3t52 [ 1 ] ; _rtB
-> o4ip5fjwmd [ 2 ] = _rtX -> chvazs3t52 [ 2 ] ; pvstwylnp5 =
muDoubleScalarCos ( _rtB -> o4ip5fjwmd [ 0 ] ) ; erv0ccpdoz =
muDoubleScalarCos ( _rtB -> o4ip5fjwmd [ 1 ] ) ; emusy1jcfj = ssGetT ( S ) ;
emusy1jcfj *= _rtP -> P_19 ; _rtB -> pu0crwpbdc = 0.4 ; _rtB -> bgwadxh2eq [
0 ] = _rtX -> p4te2bxuor [ 0 ] ; _rtB -> bgwadxh2eq [ 1 ] = _rtX ->
p4te2bxuor [ 1 ] ; _rtB -> bgwadxh2eq [ 2 ] = _rtX -> p4te2bxuor [ 2 ] ;
ieaqk3zo1a = _rtB -> pu0crwpbdc - _rtB -> bgwadxh2eq [ 2 ] ; fegfldezgg [ 0 ]
= _rtX -> daauadhpee [ 0 ] ; fegfldezgg [ 1 ] = _rtX -> daauadhpee [ 1 ] ;
fegfldezgg [ 2 ] = _rtX -> daauadhpee [ 2 ] ; _rtB -> bz1hksabsv [ 0 ] = _rtX
-> daauadhpee [ 0 ] ; _rtB -> bz1hksabsv [ 1 ] = _rtX -> daauadhpee [ 1 ] ;
_rtB -> bz1hksabsv [ 2 ] = _rtX -> daauadhpee [ 2 ] ; if ( ( _rtDW ->
pn3mnuc0r2 >= ssGetT ( S ) ) && ( _rtDW -> hjzxgmvdja >= ssGetT ( S ) ) ) {
ip43dgoh0e = 0.0 ; } else { lastTime = _rtDW -> pn3mnuc0r2 ; lastU = & _rtDW
-> kehdx0zpi5 ; if ( _rtDW -> pn3mnuc0r2 < _rtDW -> hjzxgmvdja ) { if ( _rtDW
-> hjzxgmvdja < ssGetT ( S ) ) { lastTime = _rtDW -> hjzxgmvdja ; lastU = &
_rtDW -> gqz0n4opfd ; } } else { if ( _rtDW -> pn3mnuc0r2 >= ssGetT ( S ) ) {
lastTime = _rtDW -> hjzxgmvdja ; lastU = & _rtDW -> gqz0n4opfd ; } }
ip43dgoh0e = ( _rtB -> pu0crwpbdc - * lastU ) / ( ssGetT ( S ) - lastTime ) ;
} ip43dgoh0e = ( _rtB -> bz1hksabsv [ 2 ] - ip43dgoh0e ) - _rtP -> P_22 *
ieaqk3zo1a ; _rtB -> lxbenq1gnn = ( ( ( ieaqk3zo1a + 9.81 ) - ( 10.0 *
ieaqk3zo1a + ip43dgoh0e ) * 10.0 ) - 0.3 * ip43dgoh0e ) * ( 1.65 / pvstwylnp5
/ erv0ccpdoz ) ; for ( i = 0 ; i < 4 ; i ++ ) { _rtB -> iw1wmtsjii [ i ] =
0.0 ; _rtB -> iw1wmtsjii [ i ] += _rtB -> agmy4ow3om [ i ] * _rtB ->
lxbenq1gnn ; _rtB -> iw1wmtsjii [ i ] += _rtB -> agmy4ow3om [ i + 4 ] * _rtB
-> dy0jts3w2k [ 0 ] ; _rtB -> iw1wmtsjii [ i ] += _rtB -> agmy4ow3om [ i + 8
] * _rtB -> dy0jts3w2k [ 1 ] ; _rtB -> iw1wmtsjii [ i ] += _rtB -> agmy4ow3om
[ i + 12 ] * _rtB -> dy0jts3w2k [ 2 ] ; } { real_T * * uBuffer = ( real_T * *
) & _rtDW -> ebzxnz4cux . TUbufferPtrs [ 0 ] ; real_T * * tBuffer = ( real_T
* * ) & _rtDW -> ebzxnz4cux . TUbufferPtrs [ 4 ] ; real_T simTime = ssGetT (
S ) ; real_T tMinusDelay ; { int_T i1 ; real_T * y0 = & _rtB -> kzdrzb1ky5 [
0 ] ; const real_T * u0 = & _rtB -> iw1wmtsjii [ 0 ] ; int_T * iw_Tail = &
_rtDW -> go5ejd0ffa . Tail [ 0 ] ; int_T * iw_Head = & _rtDW -> go5ejd0ffa .
Head [ 0 ] ; int_T * iw_Last = & _rtDW -> go5ejd0ffa . Last [ 0 ] ; int_T *
iw_CircularBufSize = & _rtDW -> go5ejd0ffa . CircularBufSize [ 0 ] ; for ( i1
= 0 ; i1 < 4 ; i1 ++ ) { tMinusDelay = ( ( _rtP -> P_23 > 0.0 ) ? _rtP ->
P_23 : 0.0 ) ; tMinusDelay = simTime - tMinusDelay ; if ( _rtP -> P_23 == 0.0
) y0 [ i1 ] = u0 [ i1 ] ; else y0 [ i1 ] =
Date0501_QuadModel_acc_rt_TDelayInterpolate ( tMinusDelay , 0.0 , * tBuffer ,
* uBuffer , iw_CircularBufSize [ i1 ] , & iw_Last [ i1 ] , iw_Tail [ i1 ] ,
iw_Head [ i1 ] , _rtP -> P_24 , 0 , ( boolean_T ) ( ssIsMinorTimeStep ( S )
&& ( ssGetTimeOfLastOutput ( S ) == ssGetT ( S ) ) ) ) ; tBuffer ++ ; uBuffer
++ ; } } } _rtB -> owy2fp5kxi = 0.0 ; _rtB -> owy2fp5kxi += _rtP -> P_25 *
_rtB -> kzdrzb1ky5 [ 0 ] ; _rtB -> cgewen1fqv = 0.0 ; _rtB -> cgewen1fqv +=
_rtP -> P_26 * _rtB -> kzdrzb1ky5 [ 1 ] ; _rtB -> aintmezyap = 0.0 ; _rtB ->
aintmezyap += _rtP -> P_27 * _rtB -> kzdrzb1ky5 [ 2 ] ; _rtB -> eqdle3pixo =
0.0 ; _rtB -> eqdle3pixo += _rtP -> P_28 * _rtB -> kzdrzb1ky5 [ 3 ] ;
no2x3xsfnt_idx_0 = _rtB -> owy2fp5kxi ; no2x3xsfnt_idx_1 = _rtB -> cgewen1fqv
; no2x3xsfnt_idx_2 = _rtB -> aintmezyap ; no2x3xsfnt_idx_3 = _rtB ->
eqdle3pixo ; for ( i = 0 ; i < 4 ; i ++ ) { lastTime = etztp0mvgy [ i + 12 ]
* _rtB -> eqdle3pixo + ( etztp0mvgy [ i + 8 ] * _rtB -> aintmezyap + (
etztp0mvgy [ i + 4 ] * _rtB -> cgewen1fqv + etztp0mvgy [ i ] * _rtB ->
owy2fp5kxi ) ) ; hzpktek0ad [ i ] = lastTime ; } for ( i = 0 ; i < 4 ; i ++ )
{ _rtB -> hwks2hiihf [ i ] = 0.0 ; _rtB -> hwks2hiihf [ i ] += _rtB ->
n3bbgc0hhm [ i ] * hzpktek0ad [ 0 ] ; _rtB -> hwks2hiihf [ i ] += _rtB ->
n3bbgc0hhm [ i + 4 ] * hzpktek0ad [ 1 ] ; _rtB -> hwks2hiihf [ i ] += _rtB ->
n3bbgc0hhm [ i + 8 ] * hzpktek0ad [ 2 ] ; _rtB -> hwks2hiihf [ i ] += _rtB ->
n3bbgc0hhm [ i + 12 ] * hzpktek0ad [ 3 ] ; } } if ( ssIsSampleHit ( S , 1 ,
tid ) ) { ssCallAccelRunBlock ( S , 0 , 39 , SS_CALL_MDL_OUTPUTS ) ;
ssCallAccelRunBlock ( S , 0 , 40 , SS_CALL_MDL_OUTPUTS ) ; } if (
ssIsContinuousTask ( S , tid ) ) { ip43dgoh0e = _rtP -> P_30 * _rtX ->
ez3et521ke ; _rtB -> gcwmsbqyw4 = muDoubleScalarCos ( emusy1jcfj * ip43dgoh0e
) / 4.0 ; ieaqk3zo1a = _rtB -> gcwmsbqyw4 - _rtB -> bgwadxh2eq [ 0 ] ; if ( (
_rtDW -> e2jci4idik >= ssGetT ( S ) ) && ( _rtDW -> becjzjwo4y >= ssGetT ( S
) ) ) { pvstwylnp5 = 0.0 ; } else { lastTime = _rtDW -> e2jci4idik ; lastU =
& _rtDW -> gtw34zxoo4 ; if ( _rtDW -> e2jci4idik < _rtDW -> becjzjwo4y ) { if
( _rtDW -> becjzjwo4y < ssGetT ( S ) ) { lastTime = _rtDW -> becjzjwo4y ;
lastU = & _rtDW -> chwddmtxaw ; } } else { if ( _rtDW -> e2jci4idik >= ssGetT
( S ) ) { lastTime = _rtDW -> becjzjwo4y ; lastU = & _rtDW -> chwddmtxaw ; }
} pvstwylnp5 = ( _rtB -> gcwmsbqyw4 - * lastU ) / ( ssGetT ( S ) - lastTime )
; } pvstwylnp5 = ( _rtB -> bz1hksabsv [ 0 ] - pvstwylnp5 ) - _rtP -> P_31 *
ieaqk3zo1a ; erv0ccpdoz = ( ( ieaqk3zo1a - ( pvstwylnp5 + ieaqk3zo1a ) ) -
0.1 * pvstwylnp5 ) * ( 1.65 / _rtB -> lxbenq1gnn ) ; _rtB -> mc0htwlyl2 =
muDoubleScalarSin ( ip43dgoh0e * emusy1jcfj ) / 4.0 ; ip43dgoh0e = _rtB ->
mc0htwlyl2 - _rtB -> bgwadxh2eq [ 1 ] ; if ( ( _rtDW -> g24m0yoivm >= ssGetT
( S ) ) && ( _rtDW -> o3s5v4luha >= ssGetT ( S ) ) ) { ieaqk3zo1a = 0.0 ; }
else { lastTime = _rtDW -> g24m0yoivm ; lastU = & _rtDW -> akzdemflu2 ; if (
_rtDW -> g24m0yoivm < _rtDW -> o3s5v4luha ) { if ( _rtDW -> o3s5v4luha <
ssGetT ( S ) ) { lastTime = _rtDW -> o3s5v4luha ; lastU = & _rtDW ->
ewuuifkg0s ; } } else { if ( _rtDW -> g24m0yoivm >= ssGetT ( S ) ) { lastTime
= _rtDW -> o3s5v4luha ; lastU = & _rtDW -> ewuuifkg0s ; } } ieaqk3zo1a = (
_rtB -> mc0htwlyl2 - * lastU ) / ( ssGetT ( S ) - lastTime ) ; } ieaqk3zo1a =
( _rtB -> bz1hksabsv [ 1 ] - ieaqk3zo1a ) - _rtP -> P_32 * ip43dgoh0e ;
pvstwylnp5 = ( ( ip43dgoh0e - ( ieaqk3zo1a + ip43dgoh0e ) ) - 0.1 *
ieaqk3zo1a ) * ( 1.65 / _rtB -> lxbenq1gnn ) ; ip43dgoh0e = muDoubleScalarCos
( _rtB -> o4ip5fjwmd [ 2 ] ) ; ieaqk3zo1a = muDoubleScalarSin ( _rtB ->
o4ip5fjwmd [ 2 ] ) ; u0 = ieaqk3zo1a * erv0ccpdoz - ip43dgoh0e * pvstwylnp5 ;
if ( u0 > 1.0 ) { u0 = 1.0 ; } else { if ( u0 < - 1.0 ) { u0 = - 1.0 ; } }
gmesnt5dc1 = muDoubleScalarAsin ( u0 ) ; erv0ccpdoz = ip43dgoh0e * erv0ccpdoz
+ ieaqk3zo1a * pvstwylnp5 ; u0 = erv0ccpdoz / muDoubleScalarCos ( gmesnt5dc1
) ; if ( u0 > 1.0 ) { u0 = 1.0 ; } else { if ( u0 < - 1.0 ) { u0 = - 1.0 ; }
} _rtB -> jiyfqx1lsw = muDoubleScalarAsin ( u0 ) ; emusy1jcfj =
muDoubleScalarSin ( emusy1jcfj ) ; if ( gmesnt5dc1 > _rtP -> P_33 ) { _rtB ->
owmb5rx1xj [ 0 ] = _rtP -> P_33 ; } else if ( gmesnt5dc1 < _rtP -> P_34 ) {
_rtB -> owmb5rx1xj [ 0 ] = _rtP -> P_34 ; } else { _rtB -> owmb5rx1xj [ 0 ] =
gmesnt5dc1 ; } if ( _rtB -> jiyfqx1lsw > _rtP -> P_33 ) { _rtB -> owmb5rx1xj
[ 1 ] = _rtP -> P_33 ; } else if ( _rtB -> jiyfqx1lsw < _rtP -> P_34 ) { _rtB
-> owmb5rx1xj [ 1 ] = _rtP -> P_34 ; } else { _rtB -> owmb5rx1xj [ 1 ] = _rtB
-> jiyfqx1lsw ; } if ( emusy1jcfj > _rtP -> P_33 ) { _rtB -> owmb5rx1xj [ 2 ]
= _rtP -> P_33 ; } else if ( emusy1jcfj < _rtP -> P_34 ) { _rtB -> owmb5rx1xj
[ 2 ] = _rtP -> P_34 ; } else { _rtB -> owmb5rx1xj [ 2 ] = emusy1jcfj ; } }
if ( ssIsSampleHit ( S , 1 , tid ) ) { ssCallAccelRunBlock ( S , 0 , 64 ,
SS_CALL_MDL_OUTPUTS ) ; ssCallAccelRunBlock ( S , 0 , 65 ,
SS_CALL_MDL_OUTPUTS ) ; ssCallAccelRunBlock ( S , 0 , 66 ,
SS_CALL_MDL_OUTPUTS ) ; ssCallAccelRunBlock ( S , 0 , 67 ,
SS_CALL_MDL_OUTPUTS ) ; ssCallAccelRunBlock ( S , 0 , 68 ,
SS_CALL_MDL_OUTPUTS ) ; ssCallAccelRunBlock ( S , 0 , 69 ,
SS_CALL_MDL_OUTPUTS ) ; } if ( ssIsContinuousTask ( S , tid ) ) { fegfldezgg
[ 0 ] = _rtX -> g1it3hbzj0 [ 0 ] ; fegfldezgg [ 1 ] = _rtX -> g1it3hbzj0 [ 1
] ; fegfldezgg [ 2 ] = _rtX -> g1it3hbzj0 [ 2 ] ; _rtB -> a1x1v4mypt [ 0 ] =
_rtX -> g1it3hbzj0 [ 0 ] - _rtB -> mc53goiddk [ 0 ] ; _rtB -> a1x1v4mypt [ 1
] = _rtX -> g1it3hbzj0 [ 1 ] - _rtB -> mc53goiddk [ 1 ] ; _rtB -> a1x1v4mypt
[ 2 ] = _rtX -> g1it3hbzj0 [ 2 ] - _rtB -> mc53goiddk [ 2 ] ; } if (
ssIsSampleHit ( S , 1 , tid ) ) { ssCallAccelRunBlock ( S , 0 , 73 ,
SS_CALL_MDL_OUTPUTS ) ; } if ( ssIsContinuousTask ( S , tid ) ) { _rtB ->
g1t5zjdcsl = ( ( _rtP -> P_39 * _rtB -> owmb5rx1xj [ 0 ] - _rtB -> o4ip5fjwmd
[ 0 ] ) * _rtP -> P_40 - _rtX -> lypsvf5tuh ) * _rtP -> P_42 ; _rtB ->
ke55q1g5lr = ( ( _rtP -> P_36 * _rtB -> owmb5rx1xj [ 0 ] - _rtB -> o4ip5fjwmd
[ 0 ] ) * _rtP -> P_37 + _rtX -> ituablqrnc ) + _rtB -> g1t5zjdcsl ; _rtB ->
gmsitgv5ho = ( ( _rtP -> P_46 * _rtB -> owmb5rx1xj [ 1 ] - _rtB -> o4ip5fjwmd
[ 1 ] ) * _rtP -> P_47 - _rtX -> gkvrijmbgy ) * _rtP -> P_49 ; _rtB ->
llyy4qqzrd = ( ( _rtP -> P_43 * _rtB -> owmb5rx1xj [ 1 ] - _rtB -> o4ip5fjwmd
[ 1 ] ) * _rtP -> P_44 + _rtX -> pa1uirjhny ) + _rtB -> gmsitgv5ho ; _rtB ->
jzuej111fr = ( ( _rtP -> P_53 * _rtB -> owmb5rx1xj [ 2 ] - _rtB -> o4ip5fjwmd
[ 2 ] ) * _rtP -> P_54 - _rtX -> jgpqgj3dd4 ) * _rtP -> P_56 ; _rtB ->
d2rtsqukqu = ( ( _rtP -> P_50 * _rtB -> owmb5rx1xj [ 2 ] - _rtB -> o4ip5fjwmd
[ 2 ] ) * _rtP -> P_51 + _rtX -> jjukormc3p ) + _rtB -> jzuej111fr ; } if (
ssIsSampleHit ( S , 1 , tid ) ) { ssCallAccelRunBlock ( S , 0 , 107 ,
SS_CALL_MDL_OUTPUTS ) ; ssCallAccelRunBlock ( S , 0 , 108 ,
SS_CALL_MDL_OUTPUTS ) ; } if ( ssIsSampleHit ( S , 2 , tid ) ) { memcpy ( &
diznsy5jit [ 0 ] , & _rtP -> P_57 [ 0 ] , 9U * sizeof ( real_T ) ) ; } if (
ssIsContinuousTask ( S , tid ) && ssIsSpecialSampleHit ( S , 2 , 0 , tid ) )
{ _rtB -> brmicytidr [ 0 ] = _rtB -> a1x1v4mypt [ 0 ] ; _rtB -> brmicytidr [
1 ] = _rtB -> a1x1v4mypt [ 1 ] ; _rtB -> brmicytidr [ 2 ] = _rtB ->
a1x1v4mypt [ 2 ] ; } if ( ssIsSampleHit ( S , 2 , tid ) ) { for ( i = 0 ; i <
3 ; i ++ ) { o2vmpwsi0u [ i ] = diznsy5jit [ i + 6 ] * _rtB -> brmicytidr [ 2
] + ( diznsy5jit [ i + 3 ] * _rtB -> brmicytidr [ 1 ] + diznsy5jit [ i ] *
_rtB -> brmicytidr [ 0 ] ) ; } _rtB -> nslazbbtyw [ 0 ] = _rtP -> P_58 *
o2vmpwsi0u [ 0 ] ; _rtB -> nslazbbtyw [ 1 ] = _rtP -> P_58 * o2vmpwsi0u [ 1 ]
; _rtB -> nslazbbtyw [ 2 ] = _rtP -> P_58 * o2vmpwsi0u [ 2 ] ; } if (
ssIsSampleHit ( S , 1 , tid ) ) { if ( ssIsSpecialSampleHit ( S , 2 , 1 , tid
) ) { _rtB -> jqoqfk3vdl [ 0 ] = _rtDW -> nzkmy44stj [ 0 ] ; _rtB ->
jqoqfk3vdl [ 1 ] = _rtDW -> nzkmy44stj [ 1 ] ; _rtB -> jqoqfk3vdl [ 2 ] =
_rtDW -> nzkmy44stj [ 2 ] ; } if ( _rtB -> jqoqfk3vdl [ 0 ] > _rtP -> P_60 )
{ _rtB -> djkizsiehk [ 0 ] = _rtP -> P_60 ; } else if ( _rtB -> jqoqfk3vdl [
0 ] < _rtP -> P_61 ) { _rtB -> djkizsiehk [ 0 ] = _rtP -> P_61 ; } else {
_rtB -> djkizsiehk [ 0 ] = _rtB -> jqoqfk3vdl [ 0 ] ; } if ( _rtB ->
jqoqfk3vdl [ 1 ] > _rtP -> P_60 ) { _rtB -> djkizsiehk [ 1 ] = _rtP -> P_60 ;
} else if ( _rtB -> jqoqfk3vdl [ 1 ] < _rtP -> P_61 ) { _rtB -> djkizsiehk [
1 ] = _rtP -> P_61 ; } else { _rtB -> djkizsiehk [ 1 ] = _rtB -> jqoqfk3vdl [
1 ] ; } if ( _rtB -> jqoqfk3vdl [ 2 ] > _rtP -> P_60 ) { _rtB -> djkizsiehk [
2 ] = _rtP -> P_60 ; } else if ( _rtB -> jqoqfk3vdl [ 2 ] < _rtP -> P_61 ) {
_rtB -> djkizsiehk [ 2 ] = _rtP -> P_61 ; } else { _rtB -> djkizsiehk [ 2 ] =
_rtB -> jqoqfk3vdl [ 2 ] ; } memcpy ( & _rtB -> h5ymunpfwh [ 0 ] , & _rtP ->
P_62 [ 0 ] , 9U * sizeof ( real_T ) ) ; } if ( ssIsContinuousTask ( S , tid )
) { ljqyvq454t [ 0 ] = _rtB -> ke55q1g5lr ; ljqyvq454t [ 1 ] = _rtB ->
llyy4qqzrd ; ljqyvq454t [ 2 ] = _rtB -> d2rtsqukqu ; for ( i = 0 ; i < 3 ; i
++ ) { aejymufvap [ i ] = _rtB -> h5ymunpfwh [ i + 6 ] * _rtB -> d2rtsqukqu +
( _rtB -> h5ymunpfwh [ i + 3 ] * _rtB -> llyy4qqzrd + _rtB -> h5ymunpfwh [ i
] * _rtB -> ke55q1g5lr ) ; } for ( i = 0 ; i < 3 ; i ++ ) { _rtB ->
hjl3yn2vss [ i ] = ( _rtB -> djkizsiehk [ i ] - ( ( _rtB -> h5ymunpfwh [ i +
3 ] * _rtB -> llyy4qqzrd + _rtB -> h5ymunpfwh [ i ] * _rtB -> ke55q1g5lr ) +
_rtB -> h5ymunpfwh [ i + 6 ] * _rtB -> d2rtsqukqu ) ) + gg2u4if1k1 [ i ] ; }
_rtB -> jahceq3fo1 = ( _rtB -> owmb5rx1xj [ 0 ] - _rtB -> o4ip5fjwmd [ 0 ] )
* _rtP -> P_63 ; _rtB -> euvltxwpzi = ( _rtB -> owmb5rx1xj [ 1 ] - _rtB ->
o4ip5fjwmd [ 1 ] ) * _rtP -> P_64 ; _rtB -> pc1mxrefvq = ( _rtB -> owmb5rx1xj
[ 2 ] - _rtB -> o4ip5fjwmd [ 2 ] ) * _rtP -> P_65 ; } if ( ssIsSampleHit ( S
, 1 , tid ) ) { for ( i = 0 ; i < 9 ; i ++ ) { _rtB -> o5xjbwdiaw [ i ] =
_rtP -> P_66 [ i ] ; _rtB -> cjxnnwl4ui [ i ] = _rtP -> P_67 [ i ] ; } } if (
ssIsContinuousTask ( S , tid ) ) { for ( i = 0 ; i < 3 ; i ++ ) { aejymufvap
[ i ] = _rtB -> o5xjbwdiaw [ i + 6 ] * fegfldezgg [ 2 ] + ( _rtB ->
o5xjbwdiaw [ i + 3 ] * fegfldezgg [ 1 ] + _rtB -> o5xjbwdiaw [ i ] *
fegfldezgg [ 0 ] ) ; } gg2u4if1k1 [ 0 ] += _rtB -> djkizsiehk [ 0 ] ;
gg2u4if1k1 [ 1 ] += _rtB -> djkizsiehk [ 1 ] ; lastTime = gg2u4if1k1 [ 2 ] +
_rtB -> djkizsiehk [ 2 ] ; for ( i = 0 ; i < 3 ; i ++ ) { ljqyvq454t [ i ] =
_rtB -> cjxnnwl4ui [ i + 6 ] * lastTime + ( _rtB -> cjxnnwl4ui [ i + 3 ] *
gg2u4if1k1 [ 1 ] + _rtB -> cjxnnwl4ui [ i ] * gg2u4if1k1 [ 0 ] ) ; } for ( i
= 0 ; i < 3 ; i ++ ) { tmp [ i ] = _rtB -> cjxnnwl4ui [ i + 6 ] * lastTime +
( _rtB -> cjxnnwl4ui [ i + 3 ] * gg2u4if1k1 [ 1 ] + _rtB -> cjxnnwl4ui [ i ]
* gg2u4if1k1 [ 0 ] ) ; } for ( i = 0 ; i < 3 ; i ++ ) { o2vmpwsi0u [ i ] =
_rtB -> o5xjbwdiaw [ i + 6 ] * fegfldezgg [ 2 ] + ( _rtB -> o5xjbwdiaw [ i +
3 ] * fegfldezgg [ 1 ] + _rtB -> o5xjbwdiaw [ i ] * fegfldezgg [ 0 ] ) ; }
_rtB -> dp5udq3ziy [ 0 ] = tmp [ 0 ] + o2vmpwsi0u [ 0 ] ; _rtB -> dp5udq3ziy
[ 1 ] = tmp [ 1 ] + o2vmpwsi0u [ 1 ] ; _rtB -> dp5udq3ziy [ 2 ] = tmp [ 2 ] +
o2vmpwsi0u [ 2 ] ; u0 = erv0ccpdoz / muDoubleScalarCos ( gmesnt5dc1 ) ; if (
u0 > 1.0 ) { u0 = 1.0 ; } else { if ( u0 < - 1.0 ) { u0 = - 1.0 ; } } _rtB ->
ahmdesug1u = muDoubleScalarAsin ( u0 ) ; } if ( ssIsSampleHit ( S , 1 , tid )
) { ssCallAccelRunBlock ( S , 0 , 132 , SS_CALL_MDL_OUTPUTS ) ;
ssCallAccelRunBlock ( S , 0 , 133 , SS_CALL_MDL_OUTPUTS ) ; } if (
ssIsContinuousTask ( S , tid ) ) { if ( ssGetTaskTime ( S , 0 ) < _rtP ->
P_68 ) { _rtB -> bfaifik410 = _rtP -> P_69 ; } else { _rtB -> bfaifik410 =
_rtP -> P_70 ; } no2x3xsfnt_idx_0 = _rtB -> hwks2hiihf [ 0 ] ;
no2x3xsfnt_idx_1 = _rtB -> hwks2hiihf [ 1 ] ; no2x3xsfnt_idx_2 = _rtB ->
hwks2hiihf [ 2 ] ; no2x3xsfnt_idx_3 = _rtB -> hwks2hiihf [ 3 ] ; } if (
ssIsSampleHit ( S , 1 , tid ) ) { memcpy ( & _rtB -> n0hskkclzj [ 0 ] , &
_rtP -> P_71 [ 0 ] , 9U * sizeof ( real_T ) ) ; } if ( ssIsContinuousTask ( S
, tid ) ) { for ( i = 0 ; i < 3 ; i ++ ) { aejymufvap [ i ] = _rtB ->
n0hskkclzj [ i + 6 ] * _rtB -> mc53goiddk [ 2 ] + ( _rtB -> n0hskkclzj [ i +
3 ] * _rtB -> mc53goiddk [ 1 ] + _rtB -> n0hskkclzj [ i ] * _rtB ->
mc53goiddk [ 0 ] ) ; } for ( i = 0 ; i < 3 ; i ++ ) { tmp [ i ] = _rtB ->
n0hskkclzj [ i + 6 ] * _rtB -> mc53goiddk [ 2 ] + ( _rtB -> n0hskkclzj [ i +
3 ] * _rtB -> mc53goiddk [ 1 ] + _rtB -> n0hskkclzj [ i ] * _rtB ->
mc53goiddk [ 0 ] ) ; } emusy1jcfj = _rtB -> mc53goiddk [ 2 ] * tmp [ 1 ] ;
lastTime = aejymufvap [ 0 ] ; u0 = aejymufvap [ 1 ] ; pvstwylnp5 = aejymufvap
[ 2 ] ; ieaqk3zo1a = aejymufvap [ 0 ] ; aejymufvap [ 0 ] = _rtB -> mc53goiddk
[ 1 ] * aejymufvap [ 2 ] - emusy1jcfj ; aejymufvap [ 0 ] = _rtP -> P_72 *
no2x3xsfnt_idx_1 - aejymufvap [ 0 ] ; aejymufvap [ 1 ] = _rtP -> P_73 *
no2x3xsfnt_idx_2 - ( _rtB -> mc53goiddk [ 2 ] * lastTime - _rtB -> mc53goiddk
[ 0 ] * pvstwylnp5 ) ; aejymufvap [ 2 ] = no2x3xsfnt_idx_3 - ( _rtB ->
mc53goiddk [ 0 ] * u0 - _rtB -> mc53goiddk [ 1 ] * ieaqk3zo1a ) ; } if (
ssIsSampleHit ( S , 1 , tid ) ) { memcpy ( & _rtB -> mfqwr1fj05 [ 0 ] , &
_rtP -> P_74 [ 0 ] , 9U * sizeof ( real_T ) ) ; } if ( ssIsContinuousTask ( S
, tid ) ) { for ( i = 0 ; i < 3 ; i ++ ) { ljqyvq454t [ i ] = _rtB ->
mfqwr1fj05 [ i + 6 ] * aejymufvap [ 2 ] + ( _rtB -> mfqwr1fj05 [ i + 3 ] *
aejymufvap [ 1 ] + _rtB -> mfqwr1fj05 [ i ] * aejymufvap [ 0 ] ) ; } } if (
ssIsSampleHit ( S , 1 , tid ) ) { _rtB -> mvfp5ywm1v [ 0 ] = _rtP -> P_75 [ 0
] ; _rtB -> mvfp5ywm1v [ 1 ] = _rtP -> P_75 [ 1 ] ; _rtB -> mvfp5ywm1v [ 2 ]
= _rtP -> P_75 [ 2 ] ; _rtB -> mvfp5ywm1v [ 3 ] = _rtP -> P_75 [ 3 ] ; } if (
ssIsContinuousTask ( S , tid ) ) { hzpktek0ad [ 0 ] *= _rtB -> mvfp5ywm1v [ 0
] ; hzpktek0ad [ 1 ] *= _rtB -> mvfp5ywm1v [ 1 ] ; hzpktek0ad [ 2 ] *= _rtB
-> mvfp5ywm1v [ 2 ] ; gmesnt5dc1 = ( ( hzpktek0ad [ 0 ] + hzpktek0ad [ 1 ] )
+ hzpktek0ad [ 2 ] ) + _rtB -> mvfp5ywm1v [ 3 ] * hzpktek0ad [ 3 ] ; _rtB ->
pspx2ocmic [ 0 ] = _rtP -> P_76 * gmesnt5dc1 * _rtB -> mc53goiddk [ 1 ] +
ljqyvq454t [ 0 ] ; _rtB -> pspx2ocmic [ 1 ] = _rtP -> P_77 * gmesnt5dc1 *
_rtB -> mc53goiddk [ 0 ] + ljqyvq454t [ 1 ] ; _rtB -> pspx2ocmic [ 2 ] =
ljqyvq454t [ 2 ] + 0.0 ; gmesnt5dc1 = muDoubleScalarCos ( _rtB -> o4ip5fjwmd
[ 0 ] ) ; ip43dgoh0e = muDoubleScalarSin ( _rtB -> o4ip5fjwmd [ 1 ] ) ;
ieaqk3zo1a = muDoubleScalarCos ( _rtB -> o4ip5fjwmd [ 2 ] ) ; pvstwylnp5 =
muDoubleScalarSin ( _rtB -> o4ip5fjwmd [ 0 ] ) ; erv0ccpdoz =
muDoubleScalarSin ( _rtB -> o4ip5fjwmd [ 2 ] ) ; _rtB -> dw0vmrq0nk = (
gmesnt5dc1 * ip43dgoh0e * ieaqk3zo1a + pvstwylnp5 * erv0ccpdoz ) *
no2x3xsfnt_idx_0 * _rtP -> P_78 ; _rtB -> obnmg0eymq = ( gmesnt5dc1 *
ip43dgoh0e * erv0ccpdoz - pvstwylnp5 * ieaqk3zo1a ) * no2x3xsfnt_idx_0 * _rtP
-> P_79 ; nh5etxksnt = gmesnt5dc1 * muDoubleScalarCos ( _rtB -> o4ip5fjwmd [
1 ] ) * no2x3xsfnt_idx_0 * _rtP -> P_80 ; } if ( ssIsSampleHit ( S , 1 , tid
) ) { _rtB -> h20cm2c4qg = _rtP -> P_81 ; } if ( ssIsContinuousTask ( S , tid
) ) { _rtB -> ooptlvvxs0 = _rtB -> h20cm2c4qg + nh5etxksnt ; } if (
ssIsSampleHit ( S , 1 , tid ) ) { _rtB -> mvci4gq5kr = _rtP -> P_89 * _rtDW
-> kb4sbkwk2c * _rtP -> P_94 ; _rtB -> lscuvcebct = _rtP -> P_85 * _rtDW ->
m23rkpyrxj * _rtP -> P_95 ; _rtB -> ownwjbxjpp = _rtP -> P_93 * _rtDW ->
hrcvjke4js * _rtP -> P_96 ; ssCallAccelRunBlock ( S , 0 , 191 ,
SS_CALL_MDL_OUTPUTS ) ; ssCallAccelRunBlock ( S , 0 , 192 ,
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
gqz0n4opfd ; } * lastU = _rtB -> pu0crwpbdc ; { real_T * * uBuffer = ( real_T
* * ) & _rtDW -> ebzxnz4cux . TUbufferPtrs [ 0 ] ; real_T * * tBuffer = (
real_T * * ) & _rtDW -> ebzxnz4cux . TUbufferPtrs [ 4 ] ; real_T simTime =
ssGetT ( S ) ; _rtDW -> go5ejd0ffa . Head [ 0 ] = ( ( _rtDW -> go5ejd0ffa .
Head [ 0 ] < ( _rtDW -> go5ejd0ffa . CircularBufSize [ 0 ] - 1 ) ) ? ( _rtDW
-> go5ejd0ffa . Head [ 0 ] + 1 ) : 0 ) ; if ( _rtDW -> go5ejd0ffa . Head [ 0
] == _rtDW -> go5ejd0ffa . Tail [ 0 ] ) { if ( !
Date0501_QuadModel_acc_rt_TDelayUpdateTailOrGrowBuf ( & _rtDW -> go5ejd0ffa .
CircularBufSize [ 0 ] , & _rtDW -> go5ejd0ffa . Tail [ 0 ] , & _rtDW ->
go5ejd0ffa . Head [ 0 ] , & _rtDW -> go5ejd0ffa . Last [ 0 ] , simTime - _rtP
-> P_23 , tBuffer , uBuffer , ( NULL ) , ( boolean_T ) 0 , false , & _rtDW ->
go5ejd0ffa . MaxNewBufSize ) ) { ssSetErrorStatus ( S ,
"tdelay memory allocation error" ) ; return ; } } ( * tBuffer ++ ) [ _rtDW ->
go5ejd0ffa . Head [ 0 ] ] = simTime ; ( * uBuffer ++ ) [ _rtDW -> go5ejd0ffa
. Head [ 0 ] ] = _rtB -> iw1wmtsjii [ 0 ] ; _rtDW -> go5ejd0ffa . Head [ 1 ]
= ( ( _rtDW -> go5ejd0ffa . Head [ 1 ] < ( _rtDW -> go5ejd0ffa .
CircularBufSize [ 1 ] - 1 ) ) ? ( _rtDW -> go5ejd0ffa . Head [ 1 ] + 1 ) : 0
) ; if ( _rtDW -> go5ejd0ffa . Head [ 1 ] == _rtDW -> go5ejd0ffa . Tail [ 1 ]
) { if ( ! Date0501_QuadModel_acc_rt_TDelayUpdateTailOrGrowBuf ( & _rtDW ->
go5ejd0ffa . CircularBufSize [ 1 ] , & _rtDW -> go5ejd0ffa . Tail [ 1 ] , &
_rtDW -> go5ejd0ffa . Head [ 1 ] , & _rtDW -> go5ejd0ffa . Last [ 1 ] ,
simTime - _rtP -> P_23 , tBuffer , uBuffer , ( NULL ) , ( boolean_T ) 0 ,
false , & _rtDW -> go5ejd0ffa . MaxNewBufSize ) ) { ssSetErrorStatus ( S ,
"tdelay memory allocation error" ) ; return ; } } ( * tBuffer ++ ) [ _rtDW ->
go5ejd0ffa . Head [ 1 ] ] = simTime ; ( * uBuffer ++ ) [ _rtDW -> go5ejd0ffa
. Head [ 1 ] ] = _rtB -> iw1wmtsjii [ 1 ] ; _rtDW -> go5ejd0ffa . Head [ 2 ]
= ( ( _rtDW -> go5ejd0ffa . Head [ 2 ] < ( _rtDW -> go5ejd0ffa .
CircularBufSize [ 2 ] - 1 ) ) ? ( _rtDW -> go5ejd0ffa . Head [ 2 ] + 1 ) : 0
) ; if ( _rtDW -> go5ejd0ffa . Head [ 2 ] == _rtDW -> go5ejd0ffa . Tail [ 2 ]
) { if ( ! Date0501_QuadModel_acc_rt_TDelayUpdateTailOrGrowBuf ( & _rtDW ->
go5ejd0ffa . CircularBufSize [ 2 ] , & _rtDW -> go5ejd0ffa . Tail [ 2 ] , &
_rtDW -> go5ejd0ffa . Head [ 2 ] , & _rtDW -> go5ejd0ffa . Last [ 2 ] ,
simTime - _rtP -> P_23 , tBuffer , uBuffer , ( NULL ) , ( boolean_T ) 0 ,
false , & _rtDW -> go5ejd0ffa . MaxNewBufSize ) ) { ssSetErrorStatus ( S ,
"tdelay memory allocation error" ) ; return ; } } ( * tBuffer ++ ) [ _rtDW ->
go5ejd0ffa . Head [ 2 ] ] = simTime ; ( * uBuffer ++ ) [ _rtDW -> go5ejd0ffa
. Head [ 2 ] ] = _rtB -> iw1wmtsjii [ 2 ] ; _rtDW -> go5ejd0ffa . Head [ 3 ]
= ( ( _rtDW -> go5ejd0ffa . Head [ 3 ] < ( _rtDW -> go5ejd0ffa .
CircularBufSize [ 3 ] - 1 ) ) ? ( _rtDW -> go5ejd0ffa . Head [ 3 ] + 1 ) : 0
) ; if ( _rtDW -> go5ejd0ffa . Head [ 3 ] == _rtDW -> go5ejd0ffa . Tail [ 3 ]
) { if ( ! Date0501_QuadModel_acc_rt_TDelayUpdateTailOrGrowBuf ( & _rtDW ->
go5ejd0ffa . CircularBufSize [ 3 ] , & _rtDW -> go5ejd0ffa . Tail [ 3 ] , &
_rtDW -> go5ejd0ffa . Head [ 3 ] , & _rtDW -> go5ejd0ffa . Last [ 3 ] ,
simTime - _rtP -> P_23 , tBuffer , uBuffer , ( NULL ) , ( boolean_T ) 0 ,
false , & _rtDW -> go5ejd0ffa . MaxNewBufSize ) ) { ssSetErrorStatus ( S ,
"tdelay memory allocation error" ) ; return ; } } ( * tBuffer ) [ _rtDW ->
go5ejd0ffa . Head [ 3 ] ] = simTime ; ( * uBuffer ) [ _rtDW -> go5ejd0ffa .
Head [ 3 ] ] = _rtB -> iw1wmtsjii [ 3 ] ; } } if ( ssIsContinuousTask ( S ,
tid ) ) { if ( _rtDW -> e2jci4idik == ( rtInf ) ) { _rtDW -> e2jci4idik =
ssGetT ( S ) ; lastU = & _rtDW -> gtw34zxoo4 ; } else if ( _rtDW ->
becjzjwo4y == ( rtInf ) ) { _rtDW -> becjzjwo4y = ssGetT ( S ) ; lastU = &
_rtDW -> chwddmtxaw ; } else if ( _rtDW -> e2jci4idik < _rtDW -> becjzjwo4y )
{ _rtDW -> e2jci4idik = ssGetT ( S ) ; lastU = & _rtDW -> gtw34zxoo4 ; } else
{ _rtDW -> becjzjwo4y = ssGetT ( S ) ; lastU = & _rtDW -> chwddmtxaw ; } *
lastU = _rtB -> gcwmsbqyw4 ; if ( _rtDW -> g24m0yoivm == ( rtInf ) ) { _rtDW
-> g24m0yoivm = ssGetT ( S ) ; lastU = & _rtDW -> akzdemflu2 ; } else if (
_rtDW -> o3s5v4luha == ( rtInf ) ) { _rtDW -> o3s5v4luha = ssGetT ( S ) ;
lastU = & _rtDW -> ewuuifkg0s ; } else if ( _rtDW -> g24m0yoivm < _rtDW ->
o3s5v4luha ) { _rtDW -> g24m0yoivm = ssGetT ( S ) ; lastU = & _rtDW ->
akzdemflu2 ; } else { _rtDW -> o3s5v4luha = ssGetT ( S ) ; lastU = & _rtDW ->
ewuuifkg0s ; } * lastU = _rtB -> mc0htwlyl2 ; } if ( ssIsSampleHit ( S , 1 ,
tid ) ) { sErr = GetErrorBuffer ( & _rtDW -> byjyo5m5d4 [ 0U ] ) ;
LibUpdate_Network ( & _rtDW -> byjyo5m5d4 [ 0U ] , & _rtB -> bhrzkponcd [ 0U
] , 80 ) ; if ( * sErr != 0 ) { ssSetErrorStatus ( S , sErr ) ;
ssSetStopRequested ( S , 1 ) ; } } if ( ssIsSampleHit ( S , 2 , tid ) ) {
_rtDW -> nzkmy44stj [ 0 ] = _rtB -> nslazbbtyw [ 0 ] ; _rtDW -> nzkmy44stj [
1 ] = _rtB -> nslazbbtyw [ 1 ] ; _rtDW -> nzkmy44stj [ 2 ] = _rtB ->
nslazbbtyw [ 2 ] ; } if ( ssIsSampleHit ( S , 1 , tid ) ) { _rtDW ->
m23rkpyrxj = rt_nrand_Upu32_Yd_f_pw_snf ( & _rtDW -> hugx1ib0cr ) * _rtP ->
P_83 + _rtP -> P_82 ; _rtDW -> kb4sbkwk2c = rt_nrand_Upu32_Yd_f_pw_snf ( &
_rtDW -> nqrt1fyklg ) * _rtP -> P_87 + _rtP -> P_86 ; _rtDW -> hrcvjke4js =
rt_nrand_Upu32_Yd_f_pw_snf ( & _rtDW -> btnlgyxfg3 ) * _rtP -> P_91 + _rtP ->
P_90 ; } UNUSED_PARAMETER ( tid ) ; }
#define MDL_DERIVATIVES
static void mdlDerivatives ( SimStruct * S ) { pllwawth2v * _rtB ; lfrdjjasqz
* _rtP ; isuplyng3d * _rtX ; ns1zbfbemf * _rtXdot ; _rtXdot = ( ( ns1zbfbemf
* ) ssGetdX ( S ) ) ; _rtX = ( ( isuplyng3d * ) ssGetContStates ( S ) ) ;
_rtP = ( ( lfrdjjasqz * ) ssGetDefaultParam ( S ) ) ; _rtB = ( ( pllwawth2v *
) _ssGetBlockIO ( S ) ) ; _rtXdot -> iidw152zl0 [ 0 ] = _rtB -> pspx2ocmic [
0 ] ; _rtXdot -> iidw152zl0 [ 1 ] = _rtB -> pspx2ocmic [ 1 ] ; _rtXdot ->
iidw152zl0 [ 2 ] = _rtB -> pspx2ocmic [ 2 ] ; _rtXdot -> awllp4kjam = 0.0 ;
_rtXdot -> awllp4kjam += _rtP -> P_2 * _rtX -> awllp4kjam ; _rtXdot ->
awllp4kjam += _rtB -> hjl3yn2vss [ 0 ] ; _rtXdot -> g11z1syag2 = 0.0 ;
_rtXdot -> g11z1syag2 += _rtP -> P_4 * _rtX -> g11z1syag2 ; _rtXdot ->
g11z1syag2 += _rtB -> hjl3yn2vss [ 1 ] ; _rtXdot -> bwqacranye = 0.0 ;
_rtXdot -> bwqacranye += _rtP -> P_6 * _rtX -> bwqacranye ; _rtXdot ->
bwqacranye += _rtB -> hjl3yn2vss [ 2 ] ; _rtXdot -> chvazs3t52 [ 0 ] = _rtB
-> mc53goiddk [ 0 ] ; _rtXdot -> chvazs3t52 [ 1 ] = _rtB -> mc53goiddk [ 1 ]
; _rtXdot -> chvazs3t52 [ 2 ] = _rtB -> mc53goiddk [ 2 ] ; _rtXdot ->
p4te2bxuor [ 0 ] = _rtB -> bz1hksabsv [ 0 ] ; _rtXdot -> p4te2bxuor [ 1 ] =
_rtB -> bz1hksabsv [ 1 ] ; _rtXdot -> p4te2bxuor [ 2 ] = _rtB -> bz1hksabsv [
2 ] ; _rtXdot -> daauadhpee [ 0 ] = _rtB -> dw0vmrq0nk ; _rtXdot ->
daauadhpee [ 1 ] = _rtB -> obnmg0eymq ; _rtXdot -> daauadhpee [ 2 ] = _rtB ->
ooptlvvxs0 ; _rtXdot -> ez3et521ke = 0.0 ; _rtXdot -> ez3et521ke += _rtP ->
P_29 * _rtX -> ez3et521ke ; _rtXdot -> ez3et521ke += _rtB -> bfaifik410 ;
_rtXdot -> g1it3hbzj0 [ 0 ] = _rtB -> dp5udq3ziy [ 0 ] ; _rtXdot ->
g1it3hbzj0 [ 1 ] = _rtB -> dp5udq3ziy [ 1 ] ; _rtXdot -> g1it3hbzj0 [ 2 ] =
_rtB -> dp5udq3ziy [ 2 ] ; _rtXdot -> ituablqrnc = _rtB -> jahceq3fo1 ;
_rtXdot -> lypsvf5tuh = _rtB -> g1t5zjdcsl ; _rtXdot -> pa1uirjhny = _rtB ->
euvltxwpzi ; _rtXdot -> gkvrijmbgy = _rtB -> gmsitgv5ho ; _rtXdot ->
jjukormc3p = _rtB -> pc1mxrefvq ; _rtXdot -> jgpqgj3dd4 = _rtB -> jzuej111fr
; _rtXdot -> hzmj4zwexc = 0.0 ; _rtXdot -> hzmj4zwexc += _rtP -> P_97 * _rtX
-> hzmj4zwexc ; _rtXdot -> hzmj4zwexc += _rtB -> lscuvcebct ; _rtXdot ->
gdwfbh1eyy = 0.0 ; _rtXdot -> gdwfbh1eyy += _rtP -> P_99 * _rtX -> gdwfbh1eyy
; _rtXdot -> gdwfbh1eyy += _rtB -> mvci4gq5kr ; _rtXdot -> djg03eirmm = 0.0 ;
_rtXdot -> djg03eirmm += _rtP -> P_101 * _rtX -> djg03eirmm ; _rtXdot ->
djg03eirmm += _rtB -> ownwjbxjpp ; } static void mdlInitializeSizes (
SimStruct * S ) { ssSetChecksumVal ( S , 0 , 2800197898U ) ; ssSetChecksumVal
( S , 1 , 627048364U ) ; ssSetChecksumVal ( S , 2 , 3185185295U ) ;
ssSetChecksumVal ( S , 3 , 4228458821U ) ; { mxArray * slVerStructMat = NULL
; mxArray * slStrMat = mxCreateString ( "simulink" ) ; char slVerChar [ 10 ]
; int status = mexCallMATLAB ( 1 , & slVerStructMat , 1 , & slStrMat , "ver"
) ; if ( status == 0 ) { mxArray * slVerMat = mxGetField ( slVerStructMat , 0
, "Version" ) ; if ( slVerMat == NULL ) { status = 1 ; } else { status =
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
