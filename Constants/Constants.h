#pragma once

#include <math.h>

// Math macros

const double ZERO = 0.0;
const double ONE = 1.0;
const double TWO = 2.0;
const double FOUR = 4.0;
const double HALF = 0.5;
const double QUARTER = 0.25;
const double FIFTH = 0.2;

const double PI = 4.0*atan(1.0);
const double Cn = 1.0/PI; // For SPH kernel

// Physical constants

const double GRAVITATION_CONSTANT = 6.67384e-11;
const double GAMMA = 1.41; // for hydrogen
const double GAMMA_W = 7.0; // Water

// Conversion Factors

const double METERS_TO_LYR = 1.0/9.4605284e15;
const double SECONDS_TO_YR = 1.0/31557600.0;
