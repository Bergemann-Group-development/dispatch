

MODULE unista_mod
  USE uniinv_mod
  private

  Integer, Parameter :: kdp = selected_real_kind(15)
  public :: unista
  private :: kdp
  private :: R_unista, I_unista, D_unista

  interface unista
    module procedure d_unista, r_unista, i_unista
  end interface unista

contains
  !===============================================================================
  SUBROUTINE D_unista (XDONT, NUNI)
    !   UNISTA = (Stable unique) Removes duplicates from an array,
    !            leaving unique entries in the order of their first
    !            appearance in the initial set.
    !  Michel Olagnon - Feb. 2000
    ! ----------------------------------------------------------------------------
    Real (kind=kdp), Dimension (:), Intent (InOut) :: XDONT
    Integer, Intent (Out) :: NUNI
    Integer, Dimension (Size(XDONT)) :: IWRKT
    Logical, Dimension (Size(XDONT)) :: IFMPTYT
    Integer :: ICRS
    ! ----------------------------------------------------------------------------
    Call UNIINV (XDONT, IWRKT)
    IFMPTYT = .True.
    NUNI = 0
    Do ICRS = 1, Size(XDONT)
      If (IFMPTYT(IWRKT(ICRS))) Then
        IFMPTYT(IWRKT(ICRS)) = .False.
        NUNI = NUNI + 1
        XDONT (NUNI) = XDONT (ICRS)
      End If
    End Do
    Return
  END SUBROUTINE D_unista

  !===============================================================================
  SUBROUTINE R_unista (XDONT, NUNI)
    !   UNISTA = (Stable unique) Removes duplicates from an array,
    !            leaving unique entries in the order of their first
    !            appearance in the initial set.
    !  Michel Olagnon - Feb. 2000
    ! ----------------------------------------------------------------------------
    Real, Dimension (:), Intent (InOut) :: XDONT
    Integer, Intent (Out) :: NUNI
    Integer, Dimension (Size(XDONT)) :: IWRKT
    Logical, Dimension (Size(XDONT)) :: IFMPTYT
    Integer :: ICRS
    ! --------------------------------------------------------------------------
    Call UNIINV (XDONT, IWRKT)
    IFMPTYT = .True.
    NUNI = 0
    Do ICRS = 1, Size(XDONT)
      If (IFMPTYT(IWRKT(ICRS))) Then
        IFMPTYT(IWRKT(ICRS)) = .False.
        NUNI = NUNI + 1
        XDONT (NUNI) = XDONT (ICRS)
      End If
    End Do
    Return
  END SUBROUTINE R_unista

  !===============================================================================
  SUBROUTINE I_unista (XDONT, NUNI)
    !   UNISTA = (Stable unique) Removes duplicates from an array,
    !            leaving unique entries in the order of their first
    !            appearance in the initial set.
    !  Michel Olagnon - Feb. 2000
    ! ----------------------------------------------------------------------------
    Integer, Dimension (:), Intent (InOut)  :: XDONT
    Integer, Intent (Out) :: NUNI
    Integer, Dimension (Size(XDONT)) :: IWRKT
    Logical, Dimension (Size(XDONT)) :: IFMPTYT
    Integer :: ICRS
    ! --------------------------------------------------------------------------
    Call UNIINV (XDONT, IWRKT)
    IFMPTYT = .True.
    NUNI = 0
    Do ICRS = 1, Size(XDONT)
      If (IFMPTYT(IWRKT(ICRS))) Then
        IFMPTYT(IWRKT(ICRS)) = .False.
        NUNI = NUNI + 1
        XDONT (NUNI) = XDONT (ICRS)
      End If
    End Do
    Return
  END SUBROUTINE I_unista
END MODULE unista_mod


