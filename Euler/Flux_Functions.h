#pragma once

#include "Pstate.h"
#include "Cstate.h"
#include "../Constants/Constants.h"

Cflux Roe(const Pstate& L, const Pstate& R, const Vector3D& unit_normal);
