!****m* constants/constants ===================================================!
! NAME                                                                         !
!   constants                                                                  !
!------------------------------------------------------------------------------!
! DESCRIPTION                                                                  !
!   Constants used throughout the code (mainly mathematical and scientific     !
!   constants).                                                                !
!                                                                              !
!   Module also includes functions for converting to and from atomic units.    !
!                                                                              !
!   "Natural" units are defined such that they are "natural" to the user;      !
!   not any of the standard "natural" unit systems..                           !
!                                                                              !
!   e.g:                                                                       !
!     energy      -   eV                                                       !
!     length      -   ang                                                      !
!     time        -   fs                                                       !
!     mass        -   amu                                                      !
!     temperature -   K                                                        !
!     velocity    -   ang/fs                                                   !
!     force       -   eV/ang                                                   !
!     pressure    -   GPa                                                      !
!------------------------------------------------------------------------------!
! AUTHOR                                                                       !
!   Aaron Hopkinson                                                            !
!****==========================================================================!
module constants

  implicit none

  private

  integer,          public, parameter ::  dp  = selected_real_kind(15,300)

  real(kind=dp),    public, parameter ::  pi      = 3.141592653589793238462643383279_dp
  real(kind=dp),    public, parameter ::  two_pi   = 2.0_dp * pi
  real(kind=dp),    public, parameter ::  halfpi  = 0.5_dp * pi
  complex(kind=dp), public, parameter ::  cmplx_i = (0.0_dp, 1.0_dp)

  integer,                                  public, parameter ::  nelements = 118
  character(len=3),  dimension(nelements),  public, parameter ::  element_symbol = (/ 'H  ', 'He ', 'Li ', 'Be ', 'B  ',  &
                                                            &    'C  ', 'N  ', 'O  ', 'F  ', 'Ne ', 'Na ', 'Mg ', 'Al ',  &
                                                            &    'Si ', 'P  ', 'S  ', 'Cl ', 'Ar ', 'K  ', 'Ca ', 'Sc ',  &
                                                            &    'Ti ', 'V  ', 'Cr ', 'Mn ', 'Fe ', 'Co ', 'Ni ', 'Cu ',  &
                                                            &    'Zn ', 'Ga ', 'Ge ', 'As ', 'Se ', 'Br ', 'Kr ', 'Rb ',  &
                                                            &    'Sr ', 'Y  ', 'Zr ', 'Nb ', 'Mo ', 'Tc ', 'Ru ', 'Rh ',  &
                                                            &    'Pd ', 'Ag ', 'Cd ', 'In ', 'Sn ', 'Sb ', 'Te ', 'I  ',  &
                                                            &    'Xe ', 'Cs ', 'Ba ', 'La ', 'Ce ', 'Pr ', 'Nd ', 'Pm ',  &
                                                            &    'Sm ', 'Eu ', 'Gd ', 'Tb ', 'Dy ', 'Ho ', 'Er ', 'Tm ',  &
                                                            &    'Yb ', 'Lu ', 'Hf ', 'Ta ', 'W  ', 'Re ', 'Os ', 'Ir ',  &
                                                            &    'Pt ', 'Au ', 'Hg ', 'Tl ', 'Pb ', 'Bi ', 'Po ', 'At ',  &
                                                            &    'Rn ', 'Fr ', 'Ra ', 'Ac ', 'Th ', 'Pa ', 'U  ', 'Np ',  &
                                                            &    'Pu ', 'Am ', 'Cm ', 'Bk ', 'Cf ', 'Es ', 'Fm ', 'Md ',  &
                                                            &    'No ', 'Lr ', 'Rf ', 'Db ', 'Sg ', 'Bh ', 'Hs ', 'Mt ',  &
                                                            &    'Ds ', 'Rg ', 'Cn ', 'Uut', 'Fl ', 'Uup', 'Lv ', 'Uus', 'Uuo' /)

  ! element masses taken from www.ptable.com
  real(kind=dp),  dimension(nelements), public, parameter ::  element_mass = (/ 1.00794_dp, 4.002602_dp, 6.941_dp,        &
    & 9.012182_dp, 10.811_dp, 12.0107_dp, 14.0067_dp, 15.9994_dp, 18.9984032_dp, 20.1797_dp, 22.98976928_dp, 24.3050_dp,  &
    & 16.9815386_dp, 28.0855_dp, 30.973762_dp, 32.065_dp, 35.453_dp, 39.948_dp, 39.0983_dp, 40.078_dp, 44.955912_dp,      &
    & 47.867_dp, 50.9415_dp, 51.9961_dp, 54.938045_dp, 55.845_dp, 58.933195_dp, 58.6934_dp, 63.546_dp, 65.38_dp,          &
    & 69.723_dp, 72.63_dp, 74.92160_dp, 78.96_dp, 97.904_dp, 83.798_dp, 85.4678_dp, 87.62_dp, 88.90585_dp, 91.224_dp,     &
    & 92.90638_dp, 95.96_dp, 97.9072_dp, 101.07_dp, 102.90550_dp, 106.42_dp, 107.8682_dp, 102.411_dp, 114.818_dp,         &
    & 118.710_dp, 121.760_dp, 127.60_dp, 126.90447_dp, 131.293_dp, 132.9054519_dp, 137.326_dp, 138.90547_dp, 140.116_dp,  &
    & 140.90765_dp, 144.242_dp, 145.0_dp, 150.36_dp, 151.964_dp, 157.25_dp, 158.92535_dp, 162.500_dp, 164.93032_dp,       &
    & 167.259_dp, 168.93421_dp, 173.054_dp, 174.9668_dp, 178.49_dp, 180.94788_dp, 183.84_dp, 186.207_dp, 190.23_dp,       &
    & 192.217_dp, 195.084_dp, 196.966569_dp, 200.59_dp, 204.3833_dp, 207.2_dp, 208.98040_dp, 208.9824_dp, 209.9871_dp,    &
    & 222.0176_dp, 223.0_dp, 226.0_dp, 227.0_dp, 232.03806_dp, 231.03588_dp, 238.02891_dp, 237.0_dp, 244.0_dp, 243.0_dp,  &
    & 247.0_dp, 247.0_dp, 251.0_dp, 252.0_dp, 257.0_dp, 258.0_dp, 259.0_dp, 262.0_dp, 261.0_dp, 262.0_dp, 266.0_dp,       &
    & 264.0_dp, 277.0_dp, 268.0_dp, 271.0_dp, 272.0_dp, 285.0_dp, 284.0_dp, 289.0_dp, 288.0_dp, 292.0_dp, 294.0_dp, 294.0_dp /)


  real(kind=dp),  public, parameter ::  elementary_charge   = 1.602176565E-19_dp                ! C
  real(kind=dp),  public, parameter ::  epsilon_0           = 8.854187817E-12_dp                ! F/m
  real(kind=dp),  public, parameter ::  plancks_const       = 6.62606957E-34_dp                 ! Js
  real(kind=dp),  public, parameter ::  speed_light         = 299792458.0_dp                    ! m/s
  real(kind=dp),  public, parameter ::  atomic_mass_unit    = 1.660538921E-27_dp                ! kg
  real(kind=dp),  public, parameter ::  electron_mass       = 9.10938291E-31_dp                 ! kg
  real(kind=dp),  public, parameter ::  kb                  = 1.3806488E-23_dp                  ! J/K

  real(kind=dp),  public, parameter ::  hbar                = plancks_const/two_pi              ! Js
  real(kind=dp),  public, parameter ::  coulomb_const       = 1.0_dp/(2.0_dp*two_pi*epsilon_0)  ! kg m^3 s^-2 C^-2
  real(kind=dp),  public, parameter ::  fine_struct_const   = coulomb_const*elementary_charge**2/(hbar*speed_light) ! dimensionless

  ! for unit conversions:
  real(kind=dp),  public, parameter ::  bohr_radius         = hbar/(electron_mass*speed_light*fine_struct_const)    ! in m
  real(kind=dp),  public, parameter ::  hartree_energy      = (fine_struct_const**2)*(speed_light**2)*electron_mass ! in J
  real(kind=dp),  public, parameter ::  atomic_time         = hbar/hartree_energy               ! in s
  real(kind=dp),  public, parameter ::  atomic_mass         = atomic_mass_unit/electron_mass    ! in kg
  real(kind=dp),  public, parameter ::  atomic_temperature  = hartree_energy/kb                 ! in K
  real(kind=dp),  public, parameter ::  atomic_pressure     = hartree_energy/(bohr_radius**3)   ! in Pa

  ! private variables
  real(kind=dp),  parameter         ::  femtosecond = 1.0E-15_dp  ! in s
  real(kind=dp),  parameter         ::  angstrom    = 1.0E-10_dp  ! in m
  real(kind=dp),  parameter         ::  gigapascal  = 1.0E09_dp   ! in Pa

  public  ::  element_symbol_to_z
  public  ::  units_natural_to_atomic
  public  ::  units_atomic_to_natural

contains

!****f* constants/element_symbol_to_z =========================================!
! NAME                                                                         !
!   element_symbol_to_z (PUBLIC)                                               !
!------------------------------------------------------------------------------!
! DESCRIPTION                                                                  !
!   Returns the atomic number (Z) for a given element symbol.                  !
!                                                                              !
!   -1 is returned if the element symbol does not exist in the element_symbol  !
!   array. (Therefore we must check for this after calling this function.)     !
!------------------------------------------------------------------------------!
! ARGUMENTS                                                                    !
!   character(len=*), intent(in)    ::  symbol                                 !
!==============================================================================!
function element_symbol_to_z(symbol)
  use io, only: io_str_to_ucase
  implicit none
  character(len=*), intent(in)    ::  symbol
  integer                         ::  element_symbol_to_z

  character(len=3)  ::  str       ! max length is 3
  integer           ::  ielement

  ! we dont give errors in here.. instead check in cell_from_file.. if we add errors here, can remove there
  str = trim(adjustl(symbol))
  element_symbol_to_z = -1

  ! return if we have an invalid length atomic symbol or string contains (some) invalid characters
  ! (digits being the most common error)
  if (scan(str,'0123456789') .ne. 0) then
    return
  end if

  call io_str_to_ucase(str(1:1))  ! first character always uppercase

  do ielement = 1, nelements
    if (str .eq. element_symbol(ielement)) then
      element_symbol_to_z = ielement
      exit
    end if
  end do
end function element_symbol_to_z


!****f* constants/units_natural_to_atomic =====================================!
! NAME                                                                         !
!   units_natural_to_atomic (PUBLIC)                                           !
!------------------------------------------------------------------------------!
! DESCRIPTION                                                                  !
!   Converts from the more user friendly units to atomic units.                !
!------------------------------------------------------------------------------!
! ARGUMENTS                                                                    !
!   real(kind=dp),  optional, intent(in)  ::  energy,                          !
!                                             length,                          !
!                                             mass,                            !
!                                             temperature,                     !
!                                             time,                            !
!                                             velocity,                        !
!                                             force,                           !
!                                             pressure                         !
!------------------------------------------------------------------------------!
! NOTES                                                                        !
!   Because all arguments are optional, the compiler will not check that we    !
!   are passing by keyword or that we are entering multiple arguments at once  !
!   - hence we rely on the user to call this correctly...                      !
!==============================================================================!
elemental function units_natural_to_atomic(energy, length, mass, temperature, time, velocity, force, pressure)
  implicit none
  real(kind=dp),  optional, intent(in)  ::  energy, length, mass, temperature, time, velocity, force, pressure
  real(kind=dp)                         ::  units_natural_to_atomic

  if (present(energy)) then
    ! input in eV - converts to J and then into Hartrees
    units_natural_to_atomic = energy*elementary_charge/hartree_energy

  else if (present(length)) then
    ! input in Angstroms - converts to m and then into Bohrs
    units_natural_to_atomic = length*angstrom/bohr_radius

  else if (present(mass)) then
    ! input in amu - scales such that electron_mass = 1
    units_natural_to_atomic = mass*atomic_mass

  else if (present(temperature)) then
    ! input in K - rescales such that kb = 1
    units_natural_to_atomic = temperature/atomic_temperature

  else if (present(time)) then
    ! input in fs - rescales..
    units_natural_to_atomic = time*femtosecond/atomic_time

  else if (present(velocity)) then
    ! input Angstroms/fs
    units_natural_to_atomic = velocity*angstrom/bohr_radius                   ! bohr per fs
    units_natural_to_atomic = units_natural_to_atomic*atomic_time/femtosecond ! bohr per atomic time unit

  else if (present(force)) then
    ! input eV/Angstrom - converts to hartrees/bohr
    units_natural_to_atomic = force*elementary_charge/hartree_energy          ! Hartrees per Angstrom
    units_natural_to_atomic = units_natural_to_atomic*bohr_radius/angstrom    ! Hartrees per Bohr

  else if (present(pressure)) then
    ! input in GPa
    units_natural_to_atomic = pressure*gigapascal/atomic_pressure
  end if
end function units_natural_to_atomic


!****f* constants/units_atomic_to_natural =====================================!
! NAME                                                                         !
!   units_atomic_to_natural (PUBLIC)                                           !
!------------------------------------------------------------------------------!
! DESCRIPTION                                                                  !
!   Converts from atomic units to a more user friendly set of units.           !
!------------------------------------------------------------------------------!
! ARGUMENTS                                                                    !
!   real(kind=dp),  optional, intent(in)  ::  energy,                          !
!                                             length,                          !
!                                             mass,                            !
!                                             temperature,                     !
!                                             time,                            !
!                                             velocity,                        !
!                                             force,                           !
!                                             pressure                         !
!------------------------------------------------------------------------------!
! NOTES                                                                        !
!   Because all arguments are optional, the compiler will not check that we    !
!   are passing by keyword or that we are entering multiple arguments at once  !
!   - hence we rely on the user to call this correctly...                      !
!==============================================================================!
elemental function units_atomic_to_natural(energy, length, mass, temperature, time, velocity, force, pressure)
  implicit none
  real(kind=dp),  optional, intent(in)  ::  energy, length, mass, temperature, time, velocity, force, pressure
  real(kind=dp)                         ::  units_atomic_to_natural

  if (present(energy)) then
    ! input in Hartrees - converts to J and then into eV
    units_atomic_to_natural = energy*hartree_energy/elementary_charge

  else if (present(length)) then
    ! input in Bohrs - converts to m and then into Angstroms
    units_atomic_to_natural = length*bohr_radius/angstrom

  else if (present(mass)) then
    ! rescales from electron_mass = 1 to standard amu
    units_atomic_to_natural = mass/atomic_mass

  else if (present(temperature)) then
    ! rescales such that kb is standard value - result in K
    units_atomic_to_natural = temperature*atomic_temperature

  else if (present(time)) then
    ! input in atomic time - rescales to seconds and then into fs
    units_atomic_to_natural = time*atomic_time/femtosecond

  else if (present(velocity)) then
    ! input in bohr/0.024..fs
    units_atomic_to_natural = velocity*bohr_radius/angstrom                   ! angstrom per atomic time unit
    units_atomic_to_natural = units_atomic_to_natural*femtosecond/atomic_time ! angstrom per fs

  else if (present(force)) then
    ! input is hartrees/bohr converts to eV/Angstrom
    units_atomic_to_natural = force*hartree_energy/elementary_charge          ! eV per Bohr
    units_atomic_to_natural = units_atomic_to_natural*angstrom/bohr_radius    ! eV per Angstrom

  else if (present(pressure)) then
    ! input E_h/(a_0**3), convert to Pa, then into GPa
    units_atomic_to_natural = pressure*atomic_pressure/gigapascal
  end if
end function units_atomic_to_natural
end module constants
