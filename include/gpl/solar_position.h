// Copyright (C) 2018  Brian Bailey
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, version 2.
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
#ifndef GPL_SOLAR_POSITION_H_
#define GPL_SOLAR_POSITION_H_

#include <ctime>

// Returns the Zenith angle of the sun. Assumes phi = 0 for now
double sunAngle(double latitude, double longitude, int UTC,
                struct tm time_struct);

#endif  // GPL_SOLAR_POSITION_H_
