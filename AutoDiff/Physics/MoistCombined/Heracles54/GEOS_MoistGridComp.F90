module GEOS_MoistGridComp

use qsat_util
use ras
use cloudnew
use RASPARAMS
use CLDPARAMS
use MAPL_ConstantsMod

include 'preamble'

implicit none
private
public moist_run

contains

subroutine moist_run( IM,JM,LM,DT_MOIST,IDIM,IRUN,ICMIN,DOCONVEC, &
                      CNV_PLE, PKE, PLO, PK, &
                      SIGE,IRAS,JRAS,WGT0,WGT1,CO_AUTO, &
                      RASPARAMS,CLDPARAMS,SCLMFDFR,&
                      FRLAND,KCBL,TS,KH,&
                      SEEDRAS,RASAL2_2d,CNV_FRACTION,&
                      ITRCR,irccode,CBL_TPERT,CBL_QPERT,CBL_TPERT_MXOCN,CBL_TPERT_MXLND,&
                      TH1,Q1,U1,V1,PLE,&
                      QLLS,QLCN,QILS,QICN,CLCN,CLLS)

implicit none

!Inputs - not linearized
integer, intent(in) :: IM, JM, LM, IDIM, IRUN, ICMIN, ITRCR
real(8)      , intent(in) :: DT_moist, SIGE(0:LM)
integer, intent(in), dimension(IM,JM)      :: IRAS, JRAS
real(8)      , intent(in), dimension(IM,JM,LM)   :: WGT0, WGT1
integer, intent(in), dimension(IM,JM)      :: KCBL, SEEDRAS, DOCONVEC
real(8)      , intent(in), dimension(IM,JM)      :: TS, CO_AUTO, FRLAND
real(8)      , intent(in), dimension(IM,JM,0:LM) :: PLE,KH
type(RASPARAM_TYPE), intent(in) :: RASPARAMS
type(CLDPARAM_TYPE), intent(in) :: CLDPARAMS
real(8)      , intent(in)                        :: SCLMFDFR,CBL_TPERT,CBL_QPERT,CBL_TPERT_MXOCN,CBL_TPERT_MXLND
real(8)      , intent(in), dimension(IM,JM,0:LM) :: CNV_PLE, PKE
real(8)      , intent(in), dimension(IM,JM,LM)   :: PLO, PK

!Inputs - could linearze but not be sensible, e.g. because stochasitc
real(8)      , intent(in), dimension(IM,JM)    :: RASAL2_2d,CNV_FRACTION

!Inouts
real(8)      , intent(inout), dimension(IM,JM,LM)   :: TH1,Q1,U1,V1
real(8)      , intent(inout), dimension(IM,JM,LM)   :: QLLS,QLCN,QILS,QICN,CLCN,CLLS

!Outs
integer, intent(out) :: irccode


real(8)   ,    dimension(IM,JM,LM)   :: DQS, QSS

real(8)   , dimension(IM,JM,LM)   :: QST3
real(8)   , dimension(IM,JM,LM)   :: DZET, QDDF3
real(8)   , dimension(IM,JM)      :: TPERT,QPERT

real(8)   ,    dimension(IM,JM)      :: MXDIAMx

real(8)   ,    dimension(IM,JM,LM)   :: CNV_DQLDT, CNV_MFD, CNV_PRC3, CNV_UPDF, CNV_CVW, ENTLAM

integer :: i,j,l,k


      call PRE_RASE( IM,JM,LM,TH1,Q1,TS,CNV_FRACTION,FRLAND, &
                      PLO,PKE,PK, &
                      CBL_TPERT,CBL_QPERT,CBL_TPERT_MXOCN,CBL_TPERT_MXLND, &
                      TPERT,QPERT,DQS,QSS )

      call RASE( IDIM                 , &
                 IRUN                 , &                 
                 LM                   , &
                 ICMIN                , &
                 DT_MOIST             , &
                 DOCONVEC             , &
                 MAPL8_CP              , &
                 MAPL8_ALHL            , &
                 MAPL8_ALHS            , &
                 MAPL8_TICE            , &
                 MAPL8_GRAV            , &
                 SEEDRAS              , &
                 SIGE                 , &
                 KCBL                 , &
                 WGT0                 , &
                 WGT1                 , &
                 TPERT                , &
                 QPERT                , &
                 TH1                  , &
                 Q1                   , & 
                 U1                   , &
                 V1                   , &
                 QSS                  , & 
                 DQS                  , &
                 CNV_FRACTION         , &
                 RASAL2_2d            , &
                 CO_AUTO              , &
                 CNV_PLE              , &
                 PKE                  , &
                 CNV_DQLDT            , &
                 CNV_MFD              , &
                 CNV_PRC3             , &
                 CNV_UPDF             , &
                 RASPARAMS            , &
                 ITRCR                  )


      call PRE_PROGNO_CLOUD (IM,JM,LM,TH1,PK,PLO,PKE,CNV_PLE,&
                             QST3,DZET,QDDF3,&
                             CNV_FRACTION,CLDPARAMS)


         call  PROGNO_CLOUD ( IDIM, LM         , &
                              DT_MOIST          , &
                              PLO               , &
                              CNV_PLE           , &
                              PK                , &
                              FRLAND            , &
                              KH                , &
                              CNV_MFD           , &
                              CNV_DQLDT         , &
                              CNV_PRC3          , &
                              CNV_UPDF          , &
                              U1                , &
                              V1                , & 
                              TH1               , &              
                              Q1                , &
                              QLLS              , &
                              QLCN              , &
                              QILS              , &
                              QICN              , &
                              CLCN              , &
                              CLLS              , &
                              CLDPARAMS         , &
                              SCLMFDFR          , &
                              QST3              , &
                              DZET              , &
                              QDDF3             , &
                              CNV_FRACTION      )

end subroutine moist_run

 subroutine PRE_RASE( IM,JM,LM,TH1,Q1,TS,CNV_FRACTION,FRLAND, &
                      PLO,PKE,PK, &
                      CBL_TPERT,CBL_QPERT,CBL_TPERT_MXOCN,CBL_TPERT_MXLND, &
                      TPERT,QPERT,DQS,QSS )

  implicit none

  integer, intent(in) :: IM,JM,LM
  real(8), intent(in ), dimension(IM,JM,LM)   :: TH1, Q1
  real(8), intent(in ), dimension(IM,JM)      :: TS, CNV_FRACTION, FRLAND
  real(8), intent(in )                        :: CBL_TPERT,CBL_QPERT, CBL_TPERT_MXOCN, CBL_TPERT_MXLND
  real(8), intent(in ), dimension(IM,JM,0:LM) :: PKE
  real(8), intent(in ), dimension(IM,JM,LM)   :: PLO, PK

  real(8), intent(out), dimension(IM,JM,LM)   :: DQS, QSS
  real(8), intent(out), dimension(IM,JM)      :: TPERT, QPERT

  integer                                     :: i,j,l
  real(8),              dimension(IM,JM,0:LM) :: ZLE
  real(8),              dimension(IM,JM,LM)   :: TEMP1, ZLO

   TEMP1     = TH1*PK

   ZLE(:,:,LM) = 0.
   do L=LM,1,-1
      ZLE(:,:,L-1) = TH1(:,:,L) * (1.+MAPL8_VIREPS*Q1(:,:,L))
      ZLO(:,:,L  ) = ZLE(:,:,L) + (MAPL8_CP/MAPL8_GRAV)*( PKE(:,:,L)-PK (:,:,L  ) ) * ZLE(:,:,L-1)
      ZLE(:,:,L-1) = ZLO(:,:,L) + (MAPL8_CP/MAPL8_GRAV)*( PK (:,:,L)-PKE(:,:,L-1) ) * ZLE(:,:,L-1)
   end do

   do J=1,JM
      do I=1,IM
         call DQSATpert(DQS(i,j,:),QSS(i,j,:),TEMP1(i,j,:),PLO(i,j,:),LM)
      end do
   end do

   TPERT  = ABS(CBL_TPERT) * ( TS - ( TEMP1(:,:,LM)+ MAPL8_GRAV*ZLO(:,:,LM)/MAPL8_CP )  ) 
   if (CBL_TPERT < 0) then
      ! Make TPERT 0 in areas of deep convection
      TPERT = TPERT*(1.0-CNV_FRACTION)
   endif
   QPERT  = 0.0 !CBL_QPERT * ( QSSFC - Q(:,:,LM) )      !dh: CBL_QPERT = 0.0
   TPERT  = MAX( TPERT , 0.0 )
   QPERT  = MAX( QPERT , 0.0 )

   where (FRLAND<0.1) 
      TPERT = MIN( TPERT , CBL_TPERT_MXOCN ) ! ocean
   elsewhere
      TPERT = MIN( TPERT , CBL_TPERT_MXLND ) ! land
   end where

 end subroutine PRE_RASE



!      call RASE_FAST( IDIM                 , &
!                 IRUN                 , &                 
!                 LM                   , &
!                 ICMIN                , &
!                 DT_MOIST             , &
!                 MAPL8_CP              , &
!                 MAPL8_ALHL            , &
!                 MAPL8_ALHS            , &
!                 MAPL8_TICE            , &
!                 MAPL8_GRAV            , &
!                 SEEDRAS              , &
!                 SIGE                 , &
!                 KCBL                 , &
!                 WGT0                 , &
!                 WGT1                 , &
!                 TPERT                , &
!                 QPERT                , &
!                 TH1                  , &
!                 Q1                   , & 
!                 QSS                  , & 
!                 DQS                  , &
!                 CNV_FRACTION         , &
!                 RASAL2_2d            , &
!                 CO_AUTO              , &
!                 CNV_PLE              , &
!                 PKE                  , &
!                 RASPARAMS              )
end module GEOS_MoistGridComp

