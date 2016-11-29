!
! Copyright (C) 2002-2004 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!------------------------------------------------------------------------------!
    MODULE kinds
!------------------------------------------------------------------------------!

      IMPLICIT NONE
      SAVE
! ... kind definitions
      INTEGER, PARAMETER :: DP = selected_real_kind(14,200)
      PRIVATE
      PUBLIC :: DP
!
!------------------------------------------------------------------------------!
    END MODULE kinds
!------------------------------------------------------------------------------!
