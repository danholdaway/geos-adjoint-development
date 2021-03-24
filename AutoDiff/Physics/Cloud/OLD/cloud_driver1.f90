subroutine cloud_driver ( QILS, QLLS, QICN, QLCN, QIP, QLP, Q, T)


IMPLICIT NONE

!INPUTS
real(8), intent(in) :: QILS, QLLS, QICN, QLCN, Q, T
real(8), intent(inout) :: QIP, QLP


QILSP = QIP * QILS / (QILS + QICN + 1.e-15)
QICNP = QIP * QICN / (QILS + QICN + 1.e-15)
QLLSP = QLP * QLLS / (QLLS + QLCN + 1.e-15)
QLCNP = QLP * QLCN / (QLLS + QLCN + 1.e-15)

QILSP = QILSP + 0.5 * Q + 0.7 * T
QICNP = QICNP + 0.3 * Q + 0.5 * T
QLLSP = QLLSP + 0.4 * Q + 0.6 * T
QLCNP = QLCNP + 0.2 * Q + 0.4 * T

QIP = QILSP + QICNP
QLP = QLLSP + QLCNP


end subroutine cloud_driver

