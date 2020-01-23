      ! variables for hydro Ttau relations and MLT parameters
      integer, parameter :: Nsims = 37
      integer, parameter :: Nqtau = 82
      
      character(len=2), parameter :: sNsims = '37'
      character(len=2), parameter :: sNsims1 = '38'

      real(dp), dimension(Nqtau) :: qltau, qhpf, dqhpf
      real(dp), dimension(Nsims) :: Teffsl, gravsl, FeH, alph
      real(dp), dimension(Nsims) :: sigalph, LTeff3, Logg3
      real(dp) :: aa(Nsims,Nqtau), aa1(Nsims+1,Nqtau)
      integer :: iadj(6*Nsims), iend(Nsims), isrt(Nsims)
      integer :: ist

      ! data for interpolant
      integer, parameter :: nwork=6 ! 2*pm_work_size in interp_1d_def
      ! real(dp), dimension(4*Nqtau) :: f1
      ! real(dp), dimension(Nqtau*nwork) :: work1
      real(dp), pointer :: f1(:), work1(:)
      character(len=256) :: interp_dbg_str

      real(dp) :: tau_base_hydro_Ttau, lntau_base_hydro_Ttau, logtau_base_hydro_Ttau
