! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!! This submodule provides a bitwise reproducible implementation of exp().

#include "MOM_exp.h"

submodule (MOM_intrinsic_functions) MOM_exp

use, intrinsic :: iso_fortran_env, only : int32, int64
use, intrinsic :: ieee_arithmetic, only : ieee_rint

implicit none

! Molds for transfers and numerical format queries
real, parameter :: real_mold = 0.
  !< Real mold value

integer, parameter :: int_kind &
    = merge(int64, int32, storage_size(real_mold) > storage_size(0_int32))
  !< Integer kind with the same storage size as default real
integer(kind=int_kind), parameter :: int_mold = 0
  !< Integer mold value

! Floating point layout
integer, parameter :: expbit = digits(real_mold) - 1
  !< Position of lowest exponent bit
integer, parameter :: signbit = storage_size(real_mold) - 1
  !< Position of sign bit
integer, parameter :: expwidth = signbit - expbit
  !< Number of exponent bits
integer, parameter :: expbias = maxexponent(real_mold) - 1
  !< Exponent bias

! IEEE 754 special values
integer(kind=int_kind), parameter :: pos_inf_bits &
    = ishft(2_int_kind**expwidth - 1_int_kind, expbit)
  !< IEEE +Inf bit pattern
integer(kind=int_kind), parameter :: neg_inf_bits &
    = ior(pos_inf_bits, ishft(-1_int_kind, signbit))
  !< IEEE -Inf bit pattern

! Fast integer rounding offset
real, parameter :: round_bias = 1.5 * 2_int_kind**(digits(real_mold) - 1)
  !< Binary offset used to trigger rounding of fractional values
integer(int_kind), parameter :: round_bias_bits &
    = transfer(round_bias, int_mold)
  !< Bit representation of rounding bias

contains

!> Reproducible exponential function
!!
!! Compute exp(x) with bitwise reproducibility across platforms.
module procedure exp_repro
  ! ln2 estimates
  real, parameter :: ln2 = transfer(int(z'3fe62e42fefa39ef', int_kind), real_mold)
    !< log(2): 0.693147180559945309417232... [nondim]
  real, parameter :: I_ln2 = transfer(int(z'3ff71547652b82fe', int_kind), real_mold)
    !< 1 / ln2: 1.4426950408889634073599... [nondim]

  ! Threshold for special-case handling (nonfinite, overflow, underflow).
  ! exp(x) overflows for x > ~709.78 and underflows for x < ~-745.13.
  ! We use 1024.0 as the threshold (exponent = 10); any |x| >= 1024 will
  ! definitely overflow or underflow, so we handle those specially.
  ! Values in [512, 1024) are handled by the normal path with j-bias.
  integer(int_kind), parameter :: exp_threshold = int(expbias + 10, int_kind)
    !< Exponent of 1024.0 (2^10); values with larger exponent need special handling

  ! Double-real precision of ln2/ndiv used in Cody-Waite range reduction.
  ! These are the magnitudes of glibc's negln2hiN and negln2loN constants.
  real, parameter :: ln2_ndiv_hi = transfer(int(z'3f762e42fefa0000', int_kind), real_mold)
    !< Upper precision bits of ln2 / ndiv [nondim]
  real, parameter :: ln2_ndiv_lo = transfer(int(z'3d0cf79abc9e3b3a', int_kind), real_mold)
    !< Lower precision bits of ln2 / ndiv [nondim]

  ! Table constants (may not be used)
  ! Subdivide [-1/2 ln2, 1/2 ln2] to further reduce the approximation range.
  ! This allows for a smaller polynomial at the cost of a lookup table.
  integer, parameter :: ndiv = 128
    !< Number of domain subdivisions
  real, parameter :: I_ndiv = transfer(int(z'3f80000000000000', int_kind), real_mold)
    !< 1 / ndiv [nondim]
  real, parameter :: n_ln2 = ndiv * I_ln2
    !< ndiv / ln2 [nondim]
  integer(int_kind), parameter :: index_mask = int(ndiv - 1, int_kind)
    ! TODO

  ! The following lookup tables are used to reproducibily compute 2**(j/n),
  ! where 2**(j/n) is split into real precision and a scaled residual.
  !   2**(j/ndiv) = real(2**(j/N)) * (1 + residual)

  ! real, parameter :: idiv_scale_lookup(0:ndiv-1) &
  !     = [(2.**(real(i) / real(ndiv)), i=0,ndiv-1)]
  real, parameter :: idiv_scale_lookup(0:ndiv-1) = [ &
    transfer(int(z'3ff0000000000000', int_kind), real_mold), &
    transfer(int(z'3ff0163da9fb3335', int_kind), real_mold), &
    transfer(int(z'3ff02c9a3e778061', int_kind), real_mold), &
    transfer(int(z'3ff04315e86e7f85', int_kind), real_mold), &
    transfer(int(z'3ff059b0d3158574', int_kind), real_mold), &
    transfer(int(z'3ff0706b29ddf6de', int_kind), real_mold), &
    transfer(int(z'3ff0874518759bc8', int_kind), real_mold), &
    transfer(int(z'3ff09e3ecac6f383', int_kind), real_mold), &
    transfer(int(z'3ff0b5586cf9890f', int_kind), real_mold), &
    transfer(int(z'3ff0cc922b7247f7', int_kind), real_mold), &
    transfer(int(z'3ff0e3ec32d3d1a2', int_kind), real_mold), &
    transfer(int(z'3ff0fb66affed31b', int_kind), real_mold), &
    transfer(int(z'3ff11301d0125b51', int_kind), real_mold), &
    transfer(int(z'3ff12abdc06c31cc', int_kind), real_mold), &
    transfer(int(z'3ff1429aaea92de0', int_kind), real_mold), &
    transfer(int(z'3ff15a98c8a58e51', int_kind), real_mold), &
    transfer(int(z'3ff172b83c7d517b', int_kind), real_mold), &
    transfer(int(z'3ff18af9388c8dea', int_kind), real_mold), &
    transfer(int(z'3ff1a35beb6fcb75', int_kind), real_mold), &
    transfer(int(z'3ff1bbe084045cd4', int_kind), real_mold), &
    transfer(int(z'3ff1d4873168b9aa', int_kind), real_mold), &
    transfer(int(z'3ff1ed5022fcd91d', int_kind), real_mold), &
    transfer(int(z'3ff2063b88628cd6', int_kind), real_mold), &
    transfer(int(z'3ff21f49917ddc96', int_kind), real_mold), &
    transfer(int(z'3ff2387a6e756238', int_kind), real_mold), &
    transfer(int(z'3ff251ce4fb2a63f', int_kind), real_mold), &
    transfer(int(z'3ff26b4565e27cdd', int_kind), real_mold), &
    transfer(int(z'3ff284dfe1f56381', int_kind), real_mold), &
    transfer(int(z'3ff29e9df51fdee1', int_kind), real_mold), &
    transfer(int(z'3ff2b87fd0dad990', int_kind), real_mold), &
    transfer(int(z'3ff2d285a6e4030b', int_kind), real_mold), &
    transfer(int(z'3ff2ecafa93e2f56', int_kind), real_mold), &
    transfer(int(z'3ff306fe0a31b715', int_kind), real_mold), &
    transfer(int(z'3ff32170fc4cd831', int_kind), real_mold), &
    transfer(int(z'3ff33c08b26416ff', int_kind), real_mold), &
    transfer(int(z'3ff356c55f929ff1', int_kind), real_mold), &
    transfer(int(z'3ff371a7373aa9cb', int_kind), real_mold), &
    transfer(int(z'3ff38cae6d05d866', int_kind), real_mold), &
    transfer(int(z'3ff3a7db34e59ff7', int_kind), real_mold), &
    transfer(int(z'3ff3c32dc313a8e5', int_kind), real_mold), &
    transfer(int(z'3ff3dea64c123422', int_kind), real_mold), &
    transfer(int(z'3ff3fa4504ac801c', int_kind), real_mold), &
    transfer(int(z'3ff4160a21f72e2a', int_kind), real_mold), &
    transfer(int(z'3ff431f5d950a897', int_kind), real_mold), &
    transfer(int(z'3ff44e086061892d', int_kind), real_mold), &
    transfer(int(z'3ff46a41ed1d0057', int_kind), real_mold), &
    transfer(int(z'3ff486a2b5c13cd0', int_kind), real_mold), &
    transfer(int(z'3ff4a32af0d7d3de', int_kind), real_mold), &
    transfer(int(z'3ff4bfdad5362a27', int_kind), real_mold), &
    transfer(int(z'3ff4dcb299fddd0d', int_kind), real_mold), &
    transfer(int(z'3ff4f9b2769d2ca7', int_kind), real_mold), &
    transfer(int(z'3ff516daa2cf6642', int_kind), real_mold), &
    transfer(int(z'3ff5342b569d4f82', int_kind), real_mold), &
    transfer(int(z'3ff551a4ca5d920f', int_kind), real_mold), &
    transfer(int(z'3ff56f4736b527da', int_kind), real_mold), &
    transfer(int(z'3ff58d12d497c7fd', int_kind), real_mold), &
    transfer(int(z'3ff5ab07dd485429', int_kind), real_mold), &
    transfer(int(z'3ff5c9268a5946b7', int_kind), real_mold), &
    transfer(int(z'3ff5e76f15ad2148', int_kind), real_mold), &
    transfer(int(z'3ff605e1b976dc09', int_kind), real_mold), &
    transfer(int(z'3ff6247eb03a5585', int_kind), real_mold), &
    transfer(int(z'3ff6434634ccc320', int_kind), real_mold), &
    transfer(int(z'3ff6623882552225', int_kind), real_mold), &
    transfer(int(z'3ff68155d44ca973', int_kind), real_mold), &
    transfer(int(z'3ff6a09e667f3bcd', int_kind), real_mold), &
    transfer(int(z'3ff6c012750bdabf', int_kind), real_mold), &
    transfer(int(z'3ff6dfb23c651a2f', int_kind), real_mold), &
    transfer(int(z'3ff6ff7df9519484', int_kind), real_mold), &
    transfer(int(z'3ff71f75e8ec5f74', int_kind), real_mold), &
    transfer(int(z'3ff73f9a48a58174', int_kind), real_mold), &
    transfer(int(z'3ff75feb564267c9', int_kind), real_mold), &
    transfer(int(z'3ff780694fde5d3f', int_kind), real_mold), &
    transfer(int(z'3ff7a11473eb0187', int_kind), real_mold), &
    transfer(int(z'3ff7c1ed0130c132', int_kind), real_mold), &
    transfer(int(z'3ff7e2f336cf4e62', int_kind), real_mold), &
    transfer(int(z'3ff80427543e1a12', int_kind), real_mold), &
    transfer(int(z'3ff82589994cce13', int_kind), real_mold), &
    transfer(int(z'3ff8471a4623c7ad', int_kind), real_mold), &
    transfer(int(z'3ff868d99b4492ed', int_kind), real_mold), &
    transfer(int(z'3ff88ac7d98a6699', int_kind), real_mold), &
    transfer(int(z'3ff8ace5422aa0db', int_kind), real_mold), &
    transfer(int(z'3ff8cf3216b5448c', int_kind), real_mold), &
    transfer(int(z'3ff8f1ae99157736', int_kind), real_mold), &
    transfer(int(z'3ff9145b0b91ffc6', int_kind), real_mold), &
    transfer(int(z'3ff93737b0cdc5e5', int_kind), real_mold), &
    transfer(int(z'3ff95a44cbc8520f', int_kind), real_mold), &
    transfer(int(z'3ff97d829fde4e50', int_kind), real_mold), &
    transfer(int(z'3ff9a0f170ca07ba', int_kind), real_mold), &
    transfer(int(z'3ff9c49182a3f090', int_kind), real_mold), &
    transfer(int(z'3ff9e86319e32323', int_kind), real_mold), &
    transfer(int(z'3ffa0c667b5de565', int_kind), real_mold), &
    transfer(int(z'3ffa309bec4a2d33', int_kind), real_mold), &
    transfer(int(z'3ffa5503b23e255d', int_kind), real_mold), &
    transfer(int(z'3ffa799e1330b358', int_kind), real_mold), &
    transfer(int(z'3ffa9e6b5579fdbf', int_kind), real_mold), &
    transfer(int(z'3ffac36bbfd3f37a', int_kind), real_mold), &
    transfer(int(z'3ffae89f995ad3ad', int_kind), real_mold), &
    transfer(int(z'3ffb0e07298db666', int_kind), real_mold), &
    transfer(int(z'3ffb33a2b84f15fb', int_kind), real_mold), &
    transfer(int(z'3ffb59728de5593a', int_kind), real_mold), &
    transfer(int(z'3ffb7f76f2fb5e47', int_kind), real_mold), &
    transfer(int(z'3ffba5b030a1064a', int_kind), real_mold), &
    transfer(int(z'3ffbcc1e904bc1d2', int_kind), real_mold), &
    transfer(int(z'3ffbf2c25bd71e09', int_kind), real_mold), &
    transfer(int(z'3ffc199bdd85529c', int_kind), real_mold), &
    transfer(int(z'3ffc40ab5fffd07a', int_kind), real_mold), &
    transfer(int(z'3ffc67f12e57d14b', int_kind), real_mold), &
    transfer(int(z'3ffc8f6d9406e7b5', int_kind), real_mold), &
    transfer(int(z'3ffcb720dcef9069', int_kind), real_mold), &
    transfer(int(z'3ffcdf0b555dc3fa', int_kind), real_mold), &
    transfer(int(z'3ffd072d4a07897c', int_kind), real_mold), &
    transfer(int(z'3ffd2f87080d89f2', int_kind), real_mold), &
    transfer(int(z'3ffd5818dcfba487', int_kind), real_mold), &
    transfer(int(z'3ffd80e316c98398', int_kind), real_mold), &
    transfer(int(z'3ffda9e603db3285', int_kind), real_mold), &
    transfer(int(z'3ffdd321f301b460', int_kind), real_mold), &
    transfer(int(z'3ffdfc97337b9b5f', int_kind), real_mold), &
    transfer(int(z'3ffe264614f5a129', int_kind), real_mold), &
    transfer(int(z'3ffe502ee78b3ff6', int_kind), real_mold), &
    transfer(int(z'3ffe7a51fbc74c83', int_kind), real_mold), &
    transfer(int(z'3ffea4afa2a490da', int_kind), real_mold), &
    transfer(int(z'3ffecf482d8e67f1', int_kind), real_mold), &
    transfer(int(z'3ffefa1bee615a27', int_kind), real_mold), &
    transfer(int(z'3fff252b376bba97', int_kind), real_mold), &
    transfer(int(z'3fff50765b6e4540', int_kind), real_mold), &
    transfer(int(z'3fff7bfdad9cbe14', int_kind), real_mold), &
    transfer(int(z'3fffa7c1819e90d8', int_kind), real_mold), &
    transfer(int(z'3fffd3c22b8f71f1', int_kind), real_mold) &
  ]   !< Lookup table of 2**(j/n)

  ! real, parameter :: idiv_tail_lookup(0:ndiv-1) &
  !     = [(2._real128**(real(i, real128) / real(ndiv, real128)) &
  !         / 2.**(real(i) / real(ndiv)), i=0,ndiv-1)] - 1.
  real, parameter :: idiv_tail_lookup(0:ndiv-1) = [ &
    transfer(int(z'0000000000000000', int_kind), real_mold), &
    transfer(int(z'3c9b3b4f1a88bf6e', int_kind), real_mold), &
    transfer(int(z'bc7160139cd8dc5d', int_kind), real_mold), &
    transfer(int(z'bc905e7a108766d1', int_kind), real_mold), &
    transfer(int(z'3c8cd2523567f613', int_kind), real_mold), &
    transfer(int(z'bc8bce8023f98efa', int_kind), real_mold), &
    transfer(int(z'3c60f74e61e6c861', int_kind), real_mold), &
    transfer(int(z'3c90a3e45b33d399', int_kind), real_mold), &
    transfer(int(z'3c979aa65d837b6d', int_kind), real_mold), &
    transfer(int(z'3c8eb51a92fdeffc', int_kind), real_mold), &
    transfer(int(z'3c3ebe3d702f9cd1', int_kind), real_mold), &
    transfer(int(z'bc6a033489906e0b', int_kind), real_mold), &
    transfer(int(z'bc9556522a2fbd0e', int_kind), real_mold), &
    transfer(int(z'bc5080ef8c4eea55', int_kind), real_mold), &
    transfer(int(z'bc91c923b9d5f416', int_kind), real_mold), &
    transfer(int(z'3c80d3e3e95c55af', int_kind), real_mold), &
    transfer(int(z'bc801b15eaa59348', int_kind), real_mold), &
    transfer(int(z'bc8f1ff055de323d', int_kind), real_mold), &
    transfer(int(z'3c8b898c3f1353bf', int_kind), real_mold), &
    transfer(int(z'bc96d99c7611eb26', int_kind), real_mold), &
    transfer(int(z'3c9aecf73e3a2f60', int_kind), real_mold), &
    transfer(int(z'bc8fe782cb86389d', int_kind), real_mold), &
    transfer(int(z'3c8a6f4144a6c38d', int_kind), real_mold), &
    transfer(int(z'3c807a05b0e4047d', int_kind), real_mold), &
    transfer(int(z'3c968efde3a8a894', int_kind), real_mold), &
    transfer(int(z'3c875e18f274487d', int_kind), real_mold), &
    transfer(int(z'3c80472b981fe7f2', int_kind), real_mold), &
    transfer(int(z'bc96b87b3f71085e', int_kind), real_mold), &
    transfer(int(z'3c82f7e16d09ab31', int_kind), real_mold), &
    transfer(int(z'bc3d219b1a6fbffa', int_kind), real_mold), &
    transfer(int(z'3c8b3782720c0ab4', int_kind), real_mold), &
    transfer(int(z'3c6e149289cecb8f', int_kind), real_mold), &
    transfer(int(z'3c834d754db0abb6', int_kind), real_mold), &
    transfer(int(z'3c864201e2ac744c', int_kind), real_mold), &
    transfer(int(z'3c8fdd395dd3f84a', int_kind), real_mold), &
    transfer(int(z'bc86a3803b8e5b04', int_kind), real_mold), &
    transfer(int(z'bc924aedcc4b5068', int_kind), real_mold), &
    transfer(int(z'bc9907f81b512d8e', int_kind), real_mold), &
    transfer(int(z'bc71d1e83e9436d2', int_kind), real_mold), &
    transfer(int(z'bc991919b3ce1b15', int_kind), real_mold), &
    transfer(int(z'3c859f48a72a4c6d', int_kind), real_mold), &
    transfer(int(z'bc9312607a28698a', int_kind), real_mold), &
    transfer(int(z'bc58a78f4817895b', int_kind), real_mold), &
    transfer(int(z'bc7c2c9b67499a1b', int_kind), real_mold), &
    transfer(int(z'3c4363ed60c2ac11', int_kind), real_mold), &
    transfer(int(z'3c9666093b0664ef', int_kind), real_mold), &
    transfer(int(z'3c6ecce1daa10379', int_kind), real_mold), &
    transfer(int(z'3c93ff8e3f0f1230', int_kind), real_mold), &
    transfer(int(z'3c7690cebb7aafb0', int_kind), real_mold), &
    transfer(int(z'3c931dbdeb54e077', int_kind), real_mold), &
    transfer(int(z'bc8f94340071a38e', int_kind), real_mold), &
    transfer(int(z'bc87deccdc93a349', int_kind), real_mold), &
    transfer(int(z'bc78dec6bd0f385f', int_kind), real_mold), &
    transfer(int(z'bc861246ec7b5cf6', int_kind), real_mold), &
    transfer(int(z'3c93350518fdd78e', int_kind), real_mold), &
    transfer(int(z'3c7b98b72f8a9b05', int_kind), real_mold), &
    transfer(int(z'3c9063e1e21c5409', int_kind), real_mold), &
    transfer(int(z'3c34c7855019c6ea', int_kind), real_mold), &
    transfer(int(z'3c9432e62b64c035', int_kind), real_mold), &
    transfer(int(z'bc8ce44a6199769f', int_kind), real_mold), &
    transfer(int(z'bc8c33c53bef4da8', int_kind), real_mold), &
    transfer(int(z'bc845378892be9ae', int_kind), real_mold), &
    transfer(int(z'bc93cedd78565858', int_kind), real_mold), &
    transfer(int(z'3c5710aa807e1964', int_kind), real_mold), &
    transfer(int(z'bc93b3efbf5e2228', int_kind), real_mold), &
    transfer(int(z'bc6a12ad8734b982', int_kind), real_mold), &
    transfer(int(z'bc6367efb86da9ee', int_kind), real_mold), &
    transfer(int(z'bc80dc3d54e08851', int_kind), real_mold), &
    transfer(int(z'bc781f647e5a3ecf', int_kind), real_mold), &
    transfer(int(z'bc86ee4ac08b7db0', int_kind), real_mold), &
    transfer(int(z'bc8619321e55e68a', int_kind), real_mold), &
    transfer(int(z'3c909ccb5e09d4d3', int_kind), real_mold), &
    transfer(int(z'bc7b32dcb94da51d', int_kind), real_mold), &
    transfer(int(z'3c94ecfd5467c06b', int_kind), real_mold), &
    transfer(int(z'3c65ebe1abd66c55', int_kind), real_mold), &
    transfer(int(z'bc88a1c52fb3cf42', int_kind), real_mold), &
    transfer(int(z'bc9369b6f13b3734', int_kind), real_mold), &
    transfer(int(z'bc805e843a19ff1e', int_kind), real_mold), &
    transfer(int(z'bc94d450d872576e', int_kind), real_mold), &
    transfer(int(z'3c90ad675b0e8a00', int_kind), real_mold), &
    transfer(int(z'3c8db72fc1f0eab4', int_kind), real_mold), &
    transfer(int(z'bc65b6609cc5e7ff', int_kind), real_mold), &
    transfer(int(z'3c7bf68359f35f44', int_kind), real_mold), &
    transfer(int(z'bc93091fa71e3d83', int_kind), real_mold), &
    transfer(int(z'bc5da9b88b6c1e29', int_kind), real_mold), &
    transfer(int(z'bc6c23f97c90b959', int_kind), real_mold), &
    transfer(int(z'bc92434322f4f9aa', int_kind), real_mold), &
    transfer(int(z'bc85ca6cd7668e4b', int_kind), real_mold), &
    transfer(int(z'3c71affc2b91ce27', int_kind), real_mold), &
    transfer(int(z'3c6dd235e10a73bb', int_kind), real_mold), &
    transfer(int(z'bc87c50422622263', int_kind), real_mold), &
    transfer(int(z'3c8b1c86e3e231d5', int_kind), real_mold), &
    transfer(int(z'bc91bbd1d3bcbb15', int_kind), real_mold), &
    transfer(int(z'3c90cc319cee31d2', int_kind), real_mold), &
    transfer(int(z'3c8469846e735ab3', int_kind), real_mold), &
    transfer(int(z'bc82dfcd978e9db4', int_kind), real_mold), &
    transfer(int(z'3c8c1a7792cb3387', int_kind), real_mold), &
    transfer(int(z'bc907b8f4ad1d9fa', int_kind), real_mold), &
    transfer(int(z'bc55c3d956dcaeba', int_kind), real_mold), &
    transfer(int(z'bc90a40e3da6f640', int_kind), real_mold), &
    transfer(int(z'bc68d6f438ad9334', int_kind), real_mold), &
    transfer(int(z'bc91eee26b588a35', int_kind), real_mold), &
    transfer(int(z'3c74ffd70a5fddcd', int_kind), real_mold), &
    transfer(int(z'bc91bdfbfa9298ac', int_kind), real_mold), &
    transfer(int(z'3c736eae30af0cb3', int_kind), real_mold), &
    transfer(int(z'3c8ee3325c9ffd94', int_kind), real_mold), &
    transfer(int(z'3c84e08fd10959ac', int_kind), real_mold), &
    transfer(int(z'3c63cdaf384e1a67', int_kind), real_mold), &
    transfer(int(z'3c676b2c6c921968', int_kind), real_mold), &
    transfer(int(z'bc808a1883ccb5d2', int_kind), real_mold), &
    transfer(int(z'bc8fad5d3ffffa6f', int_kind), real_mold), &
    transfer(int(z'bc900dae3875a949', int_kind), real_mold), &
    transfer(int(z'3c74a385a63d07a7', int_kind), real_mold), &
    transfer(int(z'bc82919e2040220f', int_kind), real_mold), &
    transfer(int(z'3c8e5a50d5c192ac', int_kind), real_mold), &
    transfer(int(z'3c843a59ac016b4b', int_kind), real_mold), &
    transfer(int(z'bc82d52107b43e1f', int_kind), real_mold), &
    transfer(int(z'bc892ab93b470dc9', int_kind), real_mold), &
    transfer(int(z'3c74b604603a88d3', int_kind), real_mold), &
    transfer(int(z'3c83c5ec519d7271', int_kind), real_mold), &
    transfer(int(z'bc8ff7128fd391f0', int_kind), real_mold), &
    transfer(int(z'bc8dae98e223747d', int_kind), real_mold), &
    transfer(int(z'3c8ec3bc41aa2008', int_kind), real_mold), &
    transfer(int(z'3c842b94c3a9eb32', int_kind), real_mold), &
    transfer(int(z'3c8a64a931d185ee', int_kind), real_mold), &
    transfer(int(z'bc8e37bae43be3ed', int_kind), real_mold), &
    transfer(int(z'3c77893b4d91cd9d', int_kind), real_mold), &
    transfer(int(z'3c5305c14160cc89', int_kind), real_mold) &
  ]   !< Lookup table of r, where 2^(j/N) = real64(2^(j/N)) * (1 + r)

  ! Range of K = nint(x / ln2) for which direct exponent scaling is safe.
  ! Beyond this range, a bias is applied to handle subnormals and overflow.
  ! NOTE: Fortran exponent is defined as one less than IEEE exponent.
  integer, parameter :: Kmin = minexponent(real_mold) + 1
    !< Minimum K before subnormal scaling is needed
    !! Kmin = (minexponent() - 1) + 1 (min exp(r)) (+1 safety buffer)
  integer, parameter :: Kmax = maxexponent(real_mold) - 2
    !< Maximum K before overflow scaling is needed
    !! Kmax = (maxexponent() - 1) - 0 (max exp(r)) (-1 safety buffer)
  integer(kind=int_kind), parameter :: Kbias = maxexponent(real_mold) - 2
    !< Exponent adjustment used for overflow and subnormal scaling
    !! Any bias which rescales 2**K exp(r) to O(1) works here.

  ! More table stuff

  ! Range-reduction variables
  integer :: idiv
    ! Subdivision index
  integer(kind=int_kind) :: K
    ! Nearest IEEE-rounded integer to x/ln2 [nondim]
  real :: Z
    ! Nearest integer to ndiv*K + idiv for ndiv subdivisions
    ! NOTE: Z is stored as real to avoid int-real type conversions
  integer(kind=int_kind) :: Zi
    ! Integer representation of Z
  real :: r
    ! Range-reduced input, r = x - Z ln2/N [nondim]
  real :: e
    ! Exponent of range-reduced input, e = 2**(idiv/ndiv) exp(r) [nondim]
  real :: expm1_r
    ! Approximation to exp(r) - 1 [nondim]
  real :: idiv_scale
    ! Table value for 2**(idiv/ndiv) [nondim]
  real :: idiv_tail
    ! Relative table correction for idiv_scale [nondim]

  integer(kind=int_kind) :: xb, eb
    ! Bit representations of x, e
  integer(kind=int_kind) :: x_exp
    ! Exponent bits of x

  integer(kind=int_kind) :: j
    ! Bias added to K to compensate for exponent K beyond {-1022,..,+1023}.
  integer(kind=int_kind) :: fb
    ! Bit representation 2**j, the K exponent rescale

  ! 1. Special case handling
  ! ------------------------
  ! Handle nonfinite inputs and extreme values that would overflow/underflow.
  ! Extract exponent bits and check if |x| is too large or nonfinite.

  xb = transfer(x, int_mold)
  x_exp = iand(ishft(xb, -expbit), int(2**expwidth - 1, int_kind))

  if (x_exp >= exp_threshold) then
    ! Either nonfinite (exp = 2047) or |x| >= 512
    if (x_exp >= int(2**expwidth - 1, int_kind)) then
      ! Nonfinite: +/-Inf or NaN
      ! exp(-Inf) = 0, otherwise pass-through +Inf and +/-NaN values
      ! Compute x + x to trigger `Invalid` for signaled NaNs.
      a = merge(0., x + x, xb == neg_inf_bits)
    elseif (xb >= 0) then
      ! Large positive x -> overflow to +Inf
      a = x * huge(x)
    else
      ! Large negative x -> underflow to 0
      a = tiny(x) * tiny(x)
    endif
    return
  endif

  ! 2. Range Reduction
  ! ------------------
  ! Apply a range reduction of r = x - K ln2 - (j / N) ln2, so that
  !     exp(x) = 2**K 2**(j/N) exp(r).
  ! If K = nint(x / ln2) then r is in [-1/2N ln2, 1/2N ln2] and exp(r) can be
  ! estimated by a sufficiently accurate polynomial.

  ! Compute Z = N K + j, where r = x - Z (ln 2 / N)
  Z = NEAREST_INT(x * n_ln2)
  Zi = transfer(Z + round_bias, int_mold) - round_bias_bits

  idiv = iand(Zi, index_mask)
  K = (Zi - int(idiv, int_kind)) / ndiv

  ! Since K ~ x/ln2, the terms in r = x - K ln2 will nearly cancel and there is
  ! some expected loss of precision.  To compensate, we use a Cody-Waite
  ! correction which separates ln2 into its upper 32 bits and a lower residual.

  r = (x - Z * ln2_ndiv_hi) - Z * ln2_ndiv_lo

  ! NOTE: Aggressive optimizers may reduce this to x - K * (ln2_hi + ln2_lo)
  ! which is no better than x - K ln2.

  ! 3. Polynomial approximation
  ! ---------------------------
  ! Estimate exp(r) = 1 + r * P(r) where P(r) is a 10th order Remez minimax
  ! polynomial of (exp(r) - 1)/r.  This form ensures that exp(0) = 1 exactly.

  ! Table evaluation
  idiv_scale = idiv_scale_lookup(idiv)
  idiv_tail = idiv_tail_lookup(idiv)
  expm1_r = exp_remez_expm1_estrin_4(r)
  ! Evaluate the small correction before the final addition to idiv_scale
  e = idiv_scale + idiv_scale * (idiv_tail + expm1_r)

  ! 4. Unscaling
  ! ------------
  ! Compute exp(x) = 2**K e, an exact power-of-2 calculation.
  ! Adjust scaling to compensate for subnormal output.

  ! exp(r) has range [0.707/ndiv, 1.414/ndiv], so K shifts by either 0 or -1.
  ! A resolved exponent is in the range {-1022..1023}, so K must be in
  ! {-1021,1023}.  (For safety, we reduce the range to {-1020,1022}.)

  ! Determine if K is outside the supported exponent range.
  ! If so, then apply a bias j to normalize the exponent.
  ! Kbias is chosen so that the exponent is "something near 1".
  j = merge(Kbias, 0_int_kind, K < Kmin) + merge(-Kbias, 0_int_kind, K > Kmax)

  ! Get the bit representation of exp(r)
  eb = transfer(e, int_mold)

  ! Rescale to exp(r) to exp(x), possibly including the j bias.
  eb = eb + ishft(k + j, expbit)
  a = transfer(eb, real_mold)

  ! Undo the 2**j bias as floating point multiplication.
  ! - For "normals", this has no effect.
  ! - For subnormals, this will force subnormal estimation (if enabled).
  ! - For resolvable K beyond this range, it triggers an over/underflow.
  ! - Extreme values of K have already been filtered out by the min/max step.
  fb = ishft(int(expbias, int_kind) - j, expbit)
  a = a * transfer(fb, real_mold)
end procedure exp_repro


pure function exp_remez_expm1_estrin_4(x) result(e)
  real, intent(in) :: x
    !< Input value; expected range is [-ln2/256, ln2/256]
  real :: e
    !< Approximation of exp(x) - 1 [nondim]

  real :: x2, x4
    !< Powers of x [nondim]
  real :: p01, p23, p
    !< Polynomial partial sums [nondim]

  !> glibc coefficients for exp(x) - 1 on [-ln2/256, ln2/256]
  real, parameter :: c(0:3) = [ &
      transfer(int(z'3fdffffffffffdbd', int_kind), real_mold), &
      transfer(int(z'3fc5555555555543', int_kind), real_mold), &
      transfer(int(z'3fa55555cf172b91', int_kind), real_mold), &
      transfer(int(z'3f81111167a4d017', int_kind), real_mold) &
  ]

  x2 = x * x
  x4 = x2 * x2

  p01 = c(0) + c(1) * x
  p23 = c(2) + c(3) * x

  ! Final assembly: exp(x) - 1 = x + x2*(C2+x*C3) + x4*(C4+x*C5)
  e = (x + x2 * p01) + x4 * p23
end function exp_remez_expm1_estrin_4


!> An optimized nearest-integer function for floating point reals.
!!
!! The value x is shifted from 2**K (1 + a) to 2**(digits-1+K) (1.5 + [a]),
!! causing the fractional part to be rounded according to the current IEEE
!! settings.  In almost all cases, this is nearest ties-to-even.
!!
!! The +0.5 ensures that the biased exponent of negative numbers does not drop
!! by one, which can cause half-value rounding.
!!
!! The behavior of this function does not match nint() or anint().  The nint()
!! function always ties away from zero, e.g. nint(2.5) = 3.
!!
!! It is essential that compilers not reduce (x+b)-b to x.  This can typically
!! be ensured as long as parentheses are respected.  This is managed by the
!! ENABLE_FAST_RINT macro in MOM_exp.h and assigned to NEAREST_INT().  If
!! unset, then ieee_rint() is used.
pure function fast_rint(x) result(n)
  real, intent(in) :: x
    !< Real value to be rounded to the nearest integer
  real :: n
    !< Nearest integer to x, stored as a real

  n = (x + round_bias) - round_bias
end function fast_rint

end submodule MOM_exp
