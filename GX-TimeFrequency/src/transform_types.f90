! **************************************************************************************************
!  Copyright (C) 2020-2022 Green-X library
!  This file is distributed under the terms of the APACHE2 License.
!
! **************************************************************************************************
!> \brief Defines functions used in the sine and cosine transforms
! **************************************************************************************************
module transform_types
  use kinds, only: dp

  implicit none

  private
  !> cos(iw) -> cos(it)
  integer, parameter, public :: cosine_tw = 1
  integer, parameter, public :: cosine_wt = 2
  integer, parameter, public :: sine_tw = 3

  public psi_omega, psi_tau, mat_cos_tw, mat_cos_wt, mat_sin_tw

contains

  ! psi(omega_k,x) = 2x/(x^2+omega_k^2)
  function psi_omega(omega, x_mu)
    real(dp) :: psi_omega
    real(dp), intent(in) :: omega
    real(dp), intent(in) :: x_mu

    psi_omega = 2.0_dp * x_mu / (x_mu * x_mu + omega * omega)

  end function psi_omega

  ! psi(tau_k,x) = =exp(-x*|tau_k|)
  function psi_tau(tau, x_mu)
    real(dp) :: psi_tau
    real(dp), intent(in) :: tau
    real(dp), intent(in) :: x_mu

    psi_tau = exp(-x_mu * tau)

  end function psi_tau

  ! mat_A = cos(omega_k * tau) psi(tau,x)
  function mat_cos_tw(omega, tau, x_mu)
    real(dp) :: mat_cos_tw
    real(dp), intent(in) :: omega
    real(dp), intent(in) :: tau
    real(dp), intent(in) :: x_mu

    mat_cos_tw = cos(omega * tau) * exp(-x_mu * tau)

  end function mat_cos_tw

  ! mat_A = cos(tau_k,omega) psi(omega,x)
  function mat_cos_wt(omega, tau, x_mu)
    real(dp) :: mat_cos_wt
    real(dp), intent(in) :: omega
    real(dp), intent(in) :: tau
    real(dp), intent(in) :: x_mu

    mat_cos_wt = cos(omega * tau) * psi_omega(omega, x_mu)

  end function mat_cos_wt

  ! mat_A = sin(omega_k,tau)*phi(tau,x)
  function mat_sin_tw(omega, tau, x_mu)
    real(dp) :: mat_sin_tw
    real(dp), intent(in) :: omega
    real(dp), intent(in) :: tau
    real(dp), intent(in) :: x_mu

    mat_sin_tw = sin(omega * tau) * psi_tau(tau, x_mu)

  end function mat_sin_tw

end module transform_types
