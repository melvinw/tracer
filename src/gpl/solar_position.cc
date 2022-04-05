// Copyright (C) 2018  Brian Bailey
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, version 2.
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
#include "gpl/solar_position.h"

#include <cassert>
#include <cmath>

static unsigned int julianDay(struct tm time_struct) {
  unsigned int day = time_struct.tm_mday;
  unsigned int month =
      time_struct.tm_mon + 1; // Code below assumes months from 1-12
  unsigned int year =
      time_struct.tm_year +
      1900; // Time struct contains years relative to the year 1900
  int skips_leap[] = {0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335};
  int skips_nonleap[] = {0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334};
  int *skips;

  if ((year - 2000) % 4 == 0) { // leap year
    skips = skips_leap;
  } else { // non-leap year
    skips = skips_nonleap;
  }

  return skips[month - 1] + day;
}

static double acos_safe(double x) {
  if (x < -1.0)
    x = -1.0;
  else if (x > 1.0)
    x = 1.0;
  return acosf(x);
}

static double asin_safe(double x) {
  if (x < -1.0)
    x = -1.0;
  else if (x > 1.0)
    x = 1.0;
  return asinf(x);
}

double sunAngle(double latitude, double longitude, int UTC,
                struct tm time_struct) {
  unsigned int jd = julianDay(time_struct);

  double rad = M_PI / 180.f;
  double Gamma = 2.f * M_PI * (float(jd - 1)) / 365.f;
  // solar declination angle (Iqbal Eq. 1.3.1 after Spencer)
  double delta = 0.006918f - 0.399912f * cos(Gamma) + 0.070257f * sin(Gamma) -
                 0.006758f * cos(2.f * Gamma) + 0.000907f * sin(2.f * Gamma) -
                 0.002697f * cos(3.f * Gamma) + 0.00148f * sin(3.f * Gamma);
  double EoT =
      229.18f * (0.000075f + 0.001868f * cos(Gamma) - 0.032077f * sin(Gamma) -
                 0.014615f * cos(2.f * Gamma) - 0.04089f * sin(2.f * Gamma));
  double time_dec = time_struct.tm_hour + time_struct.tm_min / 60.f; //(hours)

  int LSTM = 15.f * float(UTC); // degrees

  double TC = 4.f * (LSTM - longitude) + EoT; // minutes
  double LST = time_dec + TC / 60.f;          // hours

  double h = (LST - 12.f) * 15.f * rad; // hour angle (rad)

  // solar zenith angle
  double theta = asin_safe(sin(latitude * rad) * sin(delta) +
                           cos(latitude * rad) * cos(delta) * cos(h)); //(rad)

  assert(theta > -0.5f * M_PI && theta < 0.5f * M_PI);

  // solar elevation angle
  double phi = acos_safe((sin(delta) - sin(theta) * sin(latitude * rad)) /
                         (cos(theta) * cos(latitude * rad)));

  if (LST > 12.f) {
    phi = 2.f * M_PI - phi;
  }

  assert(phi > 0 && phi < 2.f * M_PI);

  return (90 * rad) -
         theta; // Ignoring phi for now. 90 - theta to return zenith angle
}
