/***************************************/
/* Created by FFA Flight Systems Dept. */
/* Bromma, Sweden                      */
/* 21-Dec-2000 17:39:06                */
/***************************************/

double alt_nz[3] = { 20.0000, 3000.0000, 6000.0000 };

int m_nz = 8;

double mach_nz[8] = { 0.5500, 0.8000, 0.9000, 0.9500, 1.0000, 1.0500, 1.1000, 1.2000 };

int n_nz = 3;

double nz_alpha[3][8] = {
{ 1.5000, 0.7000, 0.8000, 0.3000, 0.3000, 0.1000, 0.1400, 0.1100 },
{ 1.6000, 0.9000, 0.5500, 0.3300, 0.1400, 0.1400, 0.1000, 0.0800 },
{ 2.2000, 1.1000, 0.7000, 0.3000, 0.0500, 0.0500, 0.0500, 0.0500 } };

double nz_i[3][8] = {
{ -0.0700, -0.0318, -0.0182, -0.0350, -0.0350, -0.0300, -0.0500, -0.0500 },
{ -0.0800, -0.0410, -0.0250, -0.0385, -0.0500, -0.0500, -0.0500, -0.0500 },
{ -0.1000, -0.0501, -0.0318, -0.0400, -0.0600, -0.0600, -0.0600, -0.0600 } };

double nz_p[3][8] = {
{ -0.0045, -0.0048, -0.0027, -0.0010, -0.0010, -0.0005, -0.0005, -0.0005 },
{ -0.0090, -0.0061, -0.0037, -0.0011, -0.0005, -0.0005, -0.0005, -0.0005 },
{ -0.0150, -0.0075, -0.0048, -0.0005, -0.0002, -0.0002, -0.0002, -0.0002 } };

double nz_q[3][8] = {
{ 0.5000, 0.2418, 0.1382, 0.1200, 0.1500, 0.3000, 0.2000, 0.2000 },
{ 0.5500, 0.3110, 0.1900, 0.1320, 0.2000, 0.2000, 0.2000, 0.2000 },
{ 0.7600, 0.3801, 0.2418, 0.1700, 0.2200, 0.2200, 0.2200, 0.2200 } };

