

! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module mbd_matrix

use mbd_constants
use mbd_lapack, only: mmul, invh, invh, eigh, eigvals, eigvalsh
use mbd_utils, only: findval, exception_t, atom_index_t, is_true

use mbd_blacs, only: blacs_desc_t, blacs_all_reduce
use mbd_scalapack, only: pmmul, pinvh, pinvh, peigh, peigvalsh


use mbd_elsi, only: elsi_eigh, elsi_eigvalsh


implicit none

private
public :: contract_cross_33

type, public :: matrix_re_t
    real(dp), allocatable :: val(:, :)
    type(atom_index_t) :: idx

    type(blacs_desc_t) :: blacs

    contains
    procedure :: siz => matrix_re_siz
    procedure :: init => matrix_re_init
    procedure :: add => matrix_re_add
    procedure :: add_diag => matrix_re_add_diag
    procedure :: add_diag_scalar => matrix_re_add_diag_scalar
    procedure :: mult_cross => matrix_re_mult_cross
    procedure :: mult_rows => matrix_re_mult_rows
    procedure :: mult_cols_3n => matrix_re_mult_cols_3n
    procedure :: mult_col => matrix_re_mult_col
    procedure :: mmul => matrix_re_mmul
    procedure :: invh => matrix_re_invh
    procedure :: eigh => matrix_re_eigh
    procedure :: eigvals => matrix_re_eigvals
    procedure :: eigvalsh => matrix_re_eigvalsh
    procedure :: sum_all => matrix_re_sum_all
    procedure :: contract_n_transp => matrix_re_contract_n_transp
    procedure :: contract_n33diag_cols => matrix_re_contract_n33diag_cols
    procedure :: contract_n33_rows => matrix_re_contract_n33_rows
    procedure :: copy_from => matrix_re_copy_from
    procedure :: move_from => matrix_re_move_from
    procedure :: init_from => matrix_re_init_from
    procedure :: alloc_from => matrix_re_alloc_from
end type

type, public :: matrix_cplx_t
    complex(dp), allocatable :: val(:, :)
    type(atom_index_t) :: idx

    type(blacs_desc_t) :: blacs

    contains
    procedure :: siz => matrix_cplx_siz
    procedure :: init => matrix_cplx_init
    procedure :: add => matrix_cplx_add
    procedure :: add_diag => matrix_cplx_add_diag
    procedure :: add_diag_scalar => matrix_cplx_add_diag_scalar
    procedure :: mult_cross => matrix_cplx_mult_cross
    procedure :: mult_rows => matrix_cplx_mult_rows
    procedure :: mult_cols_3n => matrix_cplx_mult_cols_3n
    procedure :: mult_col => matrix_cplx_mult_col
    procedure :: mmul => matrix_cplx_mmul
    procedure :: eigh => matrix_cplx_eigh
    procedure :: eigvals => matrix_cplx_eigvals
    procedure :: eigvalsh => matrix_cplx_eigvalsh
    procedure :: sum_all => matrix_cplx_sum_all
    procedure :: contract_n_transp => matrix_cplx_contract_n_transp
    procedure :: contract_n33diag_cols => matrix_cplx_contract_n33diag_cols
    procedure :: contract_n33_rows => matrix_cplx_contract_n33_rows
    procedure :: copy_from => matrix_cplx_copy_from
    procedure :: move_from => matrix_cplx_move_from
    procedure :: init_from => matrix_cplx_init_from
    procedure :: alloc_from => matrix_cplx_alloc_from
end type

interface contract_cross_33
    module procedure contract_cross_33_real
    module procedure contract_cross_33_complex
end interface

contains





integer function matrix_re_siz(this, ndim) result(siz)
    class(matrix_re_t), intent(in) :: this




    integer, intent(in) :: ndim

    siz = size(this%val, ndim)
end function



subroutine matrix_re_init(this, idx, blacs)



    class(matrix_re_t), intent(out) :: this
    type(atom_index_t), intent(in) :: idx
    type(blacs_desc_t), intent(in) :: blacs

    this%idx = idx
    this%blacs = blacs
end subroutine

subroutine matrix_re_init_from(this, other)
    class(matrix_re_t), intent(out) :: this
    type(matrix_re_t), intent(in) :: other

    this%idx = other%idx
    this%blacs = other%blacs
end subroutine

subroutine matrix_re_copy_from(this, other)
    class(matrix_re_t), intent(out) :: this
    type(matrix_re_t), intent(in) :: other

    call this%init_from(other)
    this%val = other%val
end subroutine

subroutine matrix_re_move_from(this, other)
    class(matrix_re_t), intent(out) :: this
    type(matrix_re_t), intent(inout) :: other

    call this%init_from(other)
    call move_alloc(other%val, this%val)
end subroutine

subroutine matrix_re_alloc_from(this, other)
    class(matrix_re_t), intent(out) :: this
    type(matrix_re_t), intent(in) :: other

    integer :: n1, n2

    call this%init_from(other)
    n1 = other%siz(1)
    n2 = other%siz(2)
    allocate (this%val(n1, n2))
end subroutine

subroutine matrix_re_add(this, other)
    class(matrix_re_t), intent(inout) :: this
    class(matrix_re_t), intent(in) :: other

    this%val = this%val + other%val
end subroutine

subroutine matrix_re_add_diag_scalar(this, d)
    class(matrix_re_t), intent(inout) :: this
    real(dp), intent(in) :: d

    integer :: i

    call this%add_diag([(d, i = 1, this%idx%n_atoms)])
end subroutine

subroutine matrix_re_add_diag(this, d)
    class(matrix_re_t), intent(inout) :: this
    real(dp), intent(in) :: d(:)

    integer :: my_i_atom, my_j_atom, i

    do my_i_atom = 1, size(this%idx%i_atom)
        do my_j_atom = 1, size(this%idx%j_atom)
            associate ( &
                    i_atom => this%idx%i_atom(my_i_atom), &
                    j_atom => this%idx%j_atom(my_j_atom), &
                    this_diag => this%val(3*(my_i_atom-1)+1:, 3*(my_j_atom-1)+1:) &
            )
                if (i_atom /= j_atom) cycle
                do i = 1, 3
                    this_diag(i, i) = this_diag(i, i) + d(i_atom)
                end do
            end associate
        end do
    end do
end subroutine

subroutine matrix_re_mult_cross(this, b, c)
    class(matrix_re_t), intent(inout) :: this
    real(dp), intent(in) :: b(:)
    real(dp), intent(in), optional :: c(:)

    integer :: my_i_atom, my_j_atom

    do my_i_atom = 1, size(this%idx%i_atom)
        do my_j_atom = 1, size(this%idx%j_atom)
            associate ( &
                    i_atom => this%idx%i_atom(my_i_atom), &
                    j_atom => this%idx%j_atom(my_j_atom), &
                    this_sub => this%val(3*(my_i_atom-1)+1:, 3*(my_j_atom-1)+1:) &
            )
                if (present(c)) then
                    this_sub(:3, :3) = this_sub(:3, :3) * &
                        (b(i_atom)*c(j_atom)+c(i_atom)*b(j_atom))
                else
                    this_sub(:3, :3) = this_sub(:3, :3)*b(i_atom)*b(j_atom)
                end if
            end associate
        end do
    end do
end subroutine

subroutine matrix_re_mult_rows(this, b)
    class(matrix_re_t), intent(inout) :: this
    real(dp), intent(in) :: b(:)

    integer :: my_i_atom

    do my_i_atom = 1, size(this%idx%i_atom)
        associate ( &
                i_atom => this%idx%i_atom(my_i_atom), &
                this_sub => this%val(3*(my_i_atom-1)+1:, :) &
        )
            this_sub(:3, :) = this_sub(:3, :)*b(i_atom)
        end associate
    end do
end subroutine

subroutine matrix_re_mult_cols_3n(this, b)
    class(matrix_re_t), intent(inout) :: this
    real(dp), intent(in) :: b(:)

    integer :: my_j_atom, i

    do my_j_atom = 1, size(this%idx%j_atom)
        associate ( &
                b_sub => b(3*(this%idx%j_atom(my_j_atom)-1)+1:), &
                this_sub => this%val(:, 3*(my_j_atom-1)+1:) &
        )
            forall (i = 1:3) this_sub(:, i) = this_sub(:, i)*b_sub(i)
        end associate
    end do
end subroutine

subroutine matrix_re_mult_col(this, idx, a)
    class(matrix_re_t), intent(inout) :: this
    integer, intent(in) :: idx
    real(dp), intent(in) :: a(:)

    integer :: my_i_atom, my_j_atom

    do my_j_atom = 1, size(this%idx%j_atom)
        if (this%idx%j_atom(my_j_atom) /= idx) cycle
        do my_i_atom = 1, size(this%idx%i_atom)
            associate ( &
                    i_atom => this%idx%i_atom(my_i_atom), &
                    this_sub => this%val(3*(my_i_atom-1)+1:, 3*(my_j_atom-1)+1:) &
            )
                this_sub(:3, :3) = this_sub(:3, :3)*a(i_atom)
            end associate
        end do
    end do
end subroutine

subroutine matrix_re_eigh(A, eigs, exc, src, vals_only)
    class(matrix_re_t), intent(inout) :: A
    type(matrix_re_t), intent(in), optional :: src
    real(dp), intent(out) :: eigs(:)
    type(exception_t), intent(out), optional :: exc
    logical, intent(in), optional :: vals_only

    if (A%idx%parallel) then
        call elsi_eigh(A%val, A%blacs, eigs, exc, src%val, vals_only)
        return
    end if
    call eigh(A%val, eigs, exc, src%val, vals_only)
end subroutine

function matrix_re_eigvalsh(A, exc, destroy) result(eigs)
    class(matrix_re_t), intent(inout) :: A
    type(exception_t), intent(out), optional :: exc
    logical, intent(in), optional :: destroy
    real(dp) :: eigs(3*A%idx%n_atoms)

    if (A%idx%parallel) then
        eigs = elsi_eigvalsh(A%val, A%blacs, exc, destroy)
        return
    end if
    eigs = eigvalsh(A%val, exc, destroy)
end function

function matrix_re_eigvals(A, exc, destroy) result(eigs)
    class(matrix_re_t), target, intent(in) :: A
    type(exception_t), intent(out), optional :: exc
    logical, intent(in), optional :: destroy
    complex(dp) :: eigs(3*A%idx%n_atoms)

    if (A%idx%parallel) then
        exc%code = MBD_EXC_UNIMPL
        exc%msg = 'Complex general matrix diagonalization not implemented for scalapack'
    else
        eigs = eigvals(A%val, exc, destroy)
    end if
end function

real(dp) function matrix_re_sum_all(this) result(res)
    class(matrix_re_t), intent(in) :: this

    res = sum(this%val)
    if (this%idx%parallel) call blacs_all_reduce(res, this%blacs)
end function

subroutine matrix_re_contract_n_transp(this, dir, res)
    class(matrix_re_t), intent(in) :: this
    real(dp), intent(out), target :: res(:, :)
    character(len=*), intent(in) :: dir

    integer :: my_i_atom, my_j_atom
    real(dp), pointer :: res_sub(:, :)

    res(:, :) = 0d0
    do my_i_atom = 1, size(this%idx%i_atom)
        do my_j_atom = 1, size(this%idx%j_atom)
            select case (dir(1:1))
            case ('R')
                res_sub => res(:, 3*(this%idx%i_atom(my_i_atom)-1)+1:)
            case ('C')
                res_sub => res(3*(this%idx%j_atom(my_j_atom)-1)+1:, :)
            end select
            associate ( &
                    this_sub => this%val(3*(my_i_atom-1)+1:, 3*(my_j_atom-1)+1:) &
            )
                res_sub(:3, :3) = res_sub(:3, :3) + transpose(this_sub(:3, :3))
            end associate
        end do
    end do
    if (this%idx%parallel) call blacs_all_reduce(res, this%blacs)
end subroutine

function contract_cross_33_real(k_atom, A, A_prime, B, B_prime) result(res)
    type(matrix_re_t), intent(in) :: A, B
    real(dp), intent(in) :: A_prime(:, :), B_prime(:, :)
    real(dp) :: res(A%idx%n_atoms)
    integer, intent(in) :: k_atom

    integer :: my_i_atom, my_j_atom, i_atom, j_atom

    res(:) = 0d0
    my_i_atom = findval(A%idx%i_atom, k_atom)
    if (my_i_atom > 0) then
        do my_j_atom = 1, size(A%idx%j_atom)
            j_atom = A%idx%j_atom(my_j_atom)
            associate ( &
                    A_sub => A%val(3*(my_i_atom-1)+1:, 3*(my_j_atom-1)+1:), &
                    A_prime_sub => A_prime(:, 3*(j_atom-1)+1:) &
            )
                res(j_atom) = -1d0/3*sum(A_sub(:3, :3)*A_prime_sub(:, :3))
            end associate
        end do
    end if
    my_j_atom = findval(A%idx%j_atom, k_atom)
    if (my_j_atom > 0) then
        do my_i_atom = 1, size(A%idx%i_atom)
            i_atom = A%idx%i_atom(my_i_atom)
            associate ( &
                    B_sub => B%val(3*(my_i_atom-1)+1:, 3*(my_j_atom-1)+1:), &
                    B_prime_sub => B_prime(3*(i_atom-1)+1:, :) &
            )
                res(i_atom) = res(i_atom) + &
                    (-1d0/3)*sum(B_prime_sub(:3, :)*B_sub(:3, :3))
            end associate
        end do
    end if
    if (A%idx%parallel) call blacs_all_reduce(res, A%blacs)
end function

function matrix_re_contract_n33diag_cols(A) result(res)
    class(matrix_re_t), intent(in) :: A
    real(dp) :: res(A%idx%n_atoms)

    integer :: i_xyz, my_j_atom, j_atom

    res(:) = 0d0
    do my_j_atom = 1, size(A%idx%j_atom)
        j_atom = A%idx%j_atom(my_j_atom)
        do i_xyz = 1, 3
            res(j_atom) = res(j_atom) + &
                sum(A%val(i_xyz::3, 3*(my_j_atom-1)+i_xyz))
        end do
    end do
    res = res/3
    if (A%idx%parallel) call blacs_all_reduce(res, A%blacs)
end function

function matrix_re_contract_n33_rows(A) result(res)
    class(matrix_re_t), intent(in) :: A
    real(dp) :: res(A%idx%n_atoms)

    integer :: my_i_atom, i_atom

    res(:) = 0d0
    do my_i_atom = 1, size(A%idx%i_atom)
        i_atom = A%idx%i_atom(my_i_atom)
        associate (A_sub => A%val(3*(my_i_atom-1)+1:, :))
            res(i_atom) = res(i_atom) + sum(A_sub(:3, :))
        end associate
    end do
    if (A%idx%parallel) call blacs_all_reduce(res, A%blacs)
end function

type(matrix_re_t) function matrix_re_mmul( &
        A, B, transA, transB) result(C)
    class(matrix_re_t), intent(in) :: A
    type(matrix_re_t), intent(in) :: B
    character, intent(in), optional :: transA, transB

    C%idx = A%idx
    C%blacs = A%blacs
    if (.not. A%idx%parallel) then
        C%val = mmul(A%val, B%val, transA, transB)
    else
        C%val = pmmul(A%val, A%blacs, B%val, B%blacs, transA, transB, C%blacs)
    end if
end function

! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

integer function matrix_cplx_siz(this, ndim) result(siz)
    class(matrix_cplx_t), intent(in) :: this
    integer, intent(in) :: ndim

    siz = size(this%val, ndim)
end function

subroutine matrix_cplx_init(this, idx, blacs)
    class(matrix_cplx_t), intent(out) :: this
    type(atom_index_t), intent(in) :: idx
    type(blacs_desc_t), intent(in) :: blacs

    this%idx = idx
    this%blacs = blacs
end subroutine

subroutine matrix_cplx_init_from(this, other)
    class(matrix_cplx_t), intent(out) :: this
    type(matrix_cplx_t), intent(in) :: other

    this%idx = other%idx
    this%blacs = other%blacs
end subroutine

subroutine matrix_cplx_copy_from(this, other)
    class(matrix_cplx_t), intent(out) :: this
    type(matrix_cplx_t), intent(in) :: other

    call this%init_from(other)
    this%val = other%val
end subroutine

subroutine matrix_cplx_move_from(this, other)
    class(matrix_cplx_t), intent(out) :: this
    type(matrix_cplx_t), intent(inout) :: other

    call this%init_from(other)
    call move_alloc(other%val, this%val)
end subroutine

subroutine matrix_cplx_alloc_from(this, other)
    class(matrix_cplx_t), intent(out) :: this
    type(matrix_cplx_t), intent(in) :: other

    integer :: n1, n2

    call this%init_from(other)
    n1 = other%siz(1)
    n2 = other%siz(2)
    allocate (this%val(n1, n2))
end subroutine

subroutine matrix_cplx_add(this, other)
    class(matrix_cplx_t), intent(inout) :: this
    class(matrix_cplx_t), intent(in) :: other

    this%val = this%val + other%val
end subroutine

subroutine matrix_cplx_add_diag_scalar(this, d)
    class(matrix_cplx_t), intent(inout) :: this
    real(dp), intent(in) :: d

    integer :: i

    call this%add_diag([(d, i = 1, this%idx%n_atoms)])
end subroutine

subroutine matrix_cplx_add_diag(this, d)
    class(matrix_cplx_t), intent(inout) :: this
    real(dp), intent(in) :: d(:)

    integer :: my_i_atom, my_j_atom, i

    do my_i_atom = 1, size(this%idx%i_atom)
        do my_j_atom = 1, size(this%idx%j_atom)
            associate ( &
                    i_atom => this%idx%i_atom(my_i_atom), &
                    j_atom => this%idx%j_atom(my_j_atom), &
                    this_diag => this%val(3*(my_i_atom-1)+1:, 3*(my_j_atom-1)+1:) &
            )
                if (i_atom /= j_atom) cycle
                do i = 1, 3
                    this_diag(i, i) = this_diag(i, i) + d(i_atom)
                end do
            end associate
        end do
    end do
end subroutine

subroutine matrix_cplx_mult_cross(this, b, c)
    class(matrix_cplx_t), intent(inout) :: this
    real(dp), intent(in) :: b(:)
    real(dp), intent(in), optional :: c(:)

    integer :: my_i_atom, my_j_atom

    do my_i_atom = 1, size(this%idx%i_atom)
        do my_j_atom = 1, size(this%idx%j_atom)
            associate ( &
                    i_atom => this%idx%i_atom(my_i_atom), &
                    j_atom => this%idx%j_atom(my_j_atom), &
                    this_sub => this%val(3*(my_i_atom-1)+1:, 3*(my_j_atom-1)+1:) &
            )
                if (present(c)) then
                    this_sub(:3, :3) = this_sub(:3, :3) * &
                        (b(i_atom)*c(j_atom)+c(i_atom)*b(j_atom))
                else
                    this_sub(:3, :3) = this_sub(:3, :3)*b(i_atom)*b(j_atom)
                end if
            end associate
        end do
    end do
end subroutine

subroutine matrix_cplx_mult_rows(this, b)
    class(matrix_cplx_t), intent(inout) :: this
    real(dp), intent(in) :: b(:)

    integer :: my_i_atom

    do my_i_atom = 1, size(this%idx%i_atom)
        associate ( &
                i_atom => this%idx%i_atom(my_i_atom), &
                this_sub => this%val(3*(my_i_atom-1)+1:, :) &
        )
            this_sub(:3, :) = this_sub(:3, :)*b(i_atom)
        end associate
    end do
end subroutine

subroutine matrix_cplx_mult_cols_3n(this, b)
    class(matrix_cplx_t), intent(inout) :: this
    real(dp), intent(in) :: b(:)

    integer :: my_j_atom, i

    do my_j_atom = 1, size(this%idx%j_atom)
        associate ( &
                b_sub => b(3*(this%idx%j_atom(my_j_atom)-1)+1:), &
                this_sub => this%val(:, 3*(my_j_atom-1)+1:) &
        )
            forall (i = 1:3) this_sub(:, i) = this_sub(:, i)*b_sub(i)
        end associate
    end do
end subroutine

subroutine matrix_cplx_mult_col(this, idx, a)
    class(matrix_cplx_t), intent(inout) :: this
    integer, intent(in) :: idx
    real(dp), intent(in) :: a(:)

    integer :: my_i_atom, my_j_atom

    do my_j_atom = 1, size(this%idx%j_atom)
        if (this%idx%j_atom(my_j_atom) /= idx) cycle
        do my_i_atom = 1, size(this%idx%i_atom)
            associate ( &
                    i_atom => this%idx%i_atom(my_i_atom), &
                    this_sub => this%val(3*(my_i_atom-1)+1:, 3*(my_j_atom-1)+1:) &
            )
                this_sub(:3, :3) = this_sub(:3, :3)*a(i_atom)
            end associate
        end do
    end do
end subroutine

subroutine matrix_cplx_eigh(A, eigs, exc, src, vals_only)
    class(matrix_cplx_t), intent(inout) :: A
    type(matrix_cplx_t), intent(in), optional :: src
    real(dp), intent(out) :: eigs(:)
    type(exception_t), intent(out), optional :: exc
    logical, intent(in), optional :: vals_only

    if (A%idx%parallel) then
        call elsi_eigh(A%val, A%blacs, eigs, exc, src%val, vals_only)
        return
    end if
    call eigh(A%val, eigs, exc, src%val, vals_only)
end subroutine

function matrix_cplx_eigvalsh(A, exc, destroy) result(eigs)
    class(matrix_cplx_t), intent(inout) :: A
    type(exception_t), intent(out), optional :: exc
    logical, intent(in), optional :: destroy
    real(dp) :: eigs(3*A%idx%n_atoms)

    if (A%idx%parallel) then
        eigs = elsi_eigvalsh(A%val, A%blacs, exc, destroy)
        return
    end if
    eigs = eigvalsh(A%val, exc, destroy)
end function

function matrix_cplx_eigvals(A, exc, destroy) result(eigs)
    class(matrix_cplx_t), target, intent(in) :: A
    type(exception_t), intent(out), optional :: exc
    logical, intent(in), optional :: destroy
    complex(dp) :: eigs(3*A%idx%n_atoms)

    if (A%idx%parallel) then
        exc%code = MBD_EXC_UNIMPL
        exc%msg = 'Complex general matrix diagonalization not implemented for scalapack'
    else
        eigs = eigvals(A%val, exc, destroy)
    end if
end function

complex(dp) function matrix_cplx_sum_all(this) result(res)
    class(matrix_cplx_t), intent(in) :: this

    res = sum(this%val)
    if (this%idx%parallel) call blacs_all_reduce(res, this%blacs)
end function

subroutine matrix_cplx_contract_n_transp(this, dir, res)
    class(matrix_cplx_t), intent(in) :: this
    complex(dp), intent(out), target :: res(:, :)
    character(len=*), intent(in) :: dir

    integer :: my_i_atom, my_j_atom
    complex(dp), pointer :: res_sub(:, :)

    res(:, :) = 0d0
    do my_i_atom = 1, size(this%idx%i_atom)
        do my_j_atom = 1, size(this%idx%j_atom)
            select case (dir(1:1))
            case ('R')
                res_sub => res(:, 3*(this%idx%i_atom(my_i_atom)-1)+1:)
            case ('C')
                res_sub => res(3*(this%idx%j_atom(my_j_atom)-1)+1:, :)
            end select
            associate ( &
                    this_sub => this%val(3*(my_i_atom-1)+1:, 3*(my_j_atom-1)+1:) &
            )
                res_sub(:3, :3) = res_sub(:3, :3) + transpose(this_sub(:3, :3))
            end associate
        end do
    end do
    if (this%idx%parallel) call blacs_all_reduce(res, this%blacs)
end subroutine

function contract_cross_33_complex(k_atom, A, A_prime, B, B_prime) result(res)
    type(matrix_cplx_t), intent(in) :: A, B
    complex(dp), intent(in) :: A_prime(:, :), B_prime(:, :)
    complex(dp) :: res(A%idx%n_atoms)
    integer, intent(in) :: k_atom

    integer :: my_i_atom, my_j_atom, i_atom, j_atom

    res(:) = 0d0
    my_i_atom = findval(A%idx%i_atom, k_atom)
    if (my_i_atom > 0) then
        do my_j_atom = 1, size(A%idx%j_atom)
            j_atom = A%idx%j_atom(my_j_atom)
            associate ( &
                    A_sub => A%val(3*(my_i_atom-1)+1:, 3*(my_j_atom-1)+1:), &
                    A_prime_sub => A_prime(:, 3*(j_atom-1)+1:) &
            )
                res(j_atom) = -1d0/3*sum(A_sub(:3, :3)*A_prime_sub(:, :3))
            end associate
        end do
    end if
    my_j_atom = findval(A%idx%j_atom, k_atom)
    if (my_j_atom > 0) then
        do my_i_atom = 1, size(A%idx%i_atom)
            i_atom = A%idx%i_atom(my_i_atom)
            associate ( &
                    B_sub => B%val(3*(my_i_atom-1)+1:, 3*(my_j_atom-1)+1:), &
                    B_prime_sub => B_prime(3*(i_atom-1)+1:, :) &
            )
                res(i_atom) = res(i_atom) + &
                    (-1d0/3)*sum(B_prime_sub(:3, :)*B_sub(:3, :3))
            end associate
        end do
    end if
    if (A%idx%parallel) call blacs_all_reduce(res, A%blacs)
end function

function matrix_cplx_contract_n33diag_cols(A) result(res)
    class(matrix_cplx_t), intent(in) :: A
    complex(dp) :: res(A%idx%n_atoms)

    integer :: i_xyz, my_j_atom, j_atom

    res(:) = 0d0
    do my_j_atom = 1, size(A%idx%j_atom)
        j_atom = A%idx%j_atom(my_j_atom)
        do i_xyz = 1, 3
            res(j_atom) = res(j_atom) + &
                sum(A%val(i_xyz::3, 3*(my_j_atom-1)+i_xyz))
        end do
    end do
    res = res/3
    if (A%idx%parallel) call blacs_all_reduce(res, A%blacs)
end function

function matrix_cplx_contract_n33_rows(A) result(res)
    class(matrix_cplx_t), intent(in) :: A
    complex(dp) :: res(A%idx%n_atoms)

    integer :: my_i_atom, i_atom

    res(:) = 0d0
    do my_i_atom = 1, size(A%idx%i_atom)
        i_atom = A%idx%i_atom(my_i_atom)
        associate (A_sub => A%val(3*(my_i_atom-1)+1:, :))
            res(i_atom) = res(i_atom) + sum(A_sub(:3, :))
        end associate
    end do
    if (A%idx%parallel) call blacs_all_reduce(res, A%blacs)
end function

type(matrix_cplx_t) function matrix_cplx_mmul( &
        A, B, transA, transB) result(C)
    class(matrix_cplx_t), intent(in) :: A
    type(matrix_cplx_t), intent(in) :: B
    character, intent(in), optional :: transA, transB

    C%idx = A%idx
    C%blacs = A%blacs
    if (.not. A%idx%parallel) then
        C%val = mmul(A%val, B%val, transA, transB)
    else
        C%val = pmmul(A%val, A%blacs, B%val, B%blacs, transA, transB, C%blacs)
    end if
end function


subroutine matrix_re_invh(A, exc, src)
    class(matrix_re_t), intent(inout) :: A
    type(matrix_re_t), intent(in), optional :: src
    type(exception_t), intent(out), optional :: exc

    if (.not. A%idx%parallel) then
        if (present(src)) then
            call invh(A%val, exc, src%val)
        else
            call invh(A%val, exc)
        end if
    else
        if (present(src)) then
            call pinvh(A%val, A%blacs, exc, src%val)
        else
            call pinvh(A%val, A%blacs, exc)
        end if
    end if
end subroutine

end module

