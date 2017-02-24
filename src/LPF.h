// Configuration file for MCMC LPF Micrometeorite code
// Written by writeConfig.m at 2016-05-02 10:29:21

// Spacecraft Geometry
#define SC_H 8.315000e-01 // Height of spacecraft [m]
#define SC_BOT_CORNER_1_X -9.260000e-01 // x coordinate of spacecraft bottom deck corner 1 [m]
#define SC_BOT_CORNER_1_Y -2.168000e-01 // y coordinate of spacecraft bottom deck corner 1 [m]
#define SC_BOT_CORNER_2_X -9.260000e-01 // x coordinate of spacecraft bottom deck corner 2 [m]
#define SC_BOT_CORNER_2_Y 2.048000e-01 // y coordinate of spacecraft bottom deck corner 2 [m]
#define SC_BOT_CORNER_3_X -5.263000e-01 // x coordinate of spacecraft bottom deck corner 3 [m]
#define SC_BOT_CORNER_3_Y 8.970000e-01 // y coordinate of spacecraft bottom deck corner 3 [m]
#define SC_BOT_CORNER_4_X 5.163000e-01 // x coordinate of spacecraft bottom deck corner 4 [m]
#define SC_BOT_CORNER_4_Y 8.970000e-01 // y coordinate of spacecraft bottom deck corner 4 [m]
#define SC_BOT_CORNER_5_X 9.160000e-01 // x coordinate of spacecraft bottom deck corner 5 [m]
#define SC_BOT_CORNER_5_Y 2.048000e-01 // y coordinate of spacecraft bottom deck corner 5 [m]
#define SC_BOT_CORNER_6_X 9.160000e-01 // x coordinate of spacecraft bottom deck corner 6 [m]
#define SC_BOT_CORNER_6_Y -2.168000e-01 // y coordinate of spacecraft bottom deck corner 6 [m]
#define SC_BOT_CORNER_7_X 5.163000e-01 // x coordinate of spacecraft bottom deck corner 7 [m]
#define SC_BOT_CORNER_7_Y -9.090000e-01 // y coordinate of spacecraft bottom deck corner 7 [m]
#define SC_BOT_CORNER_8_X -5.263000e-01 // x coordinate of spacecraft bottom deck corner 8 [m]
#define SC_BOT_CORNER_8_Y -9.090000e-01 // y coordinate of spacecraft bottom deck corner 8 [m]
#define EOM_RB_X 1.303672e-03 // X of S/C CoM in M Frame => this defines the B frame [m]
#define EOM_RB_Y 2.370536e-03 // Y of S/C CoM in M Frame => this defines the B frame [m]
#define EOM_RB_Z 4.913037e-01 // Z of S/C CoM in M Frame => this defines the B frame [m]

// Housing 1 Geometry
#define EOM_H1SC_X 1.880000e-01 // X of housing 1 [m]
#define EOM_H1SC_Y 0.000000e+00 // Y of housing 1 [m]
#define EOM_H1SC_Z 5.325000e-01 // Z of housing 1 [m]

// Housing 2 Geometry
#define EOM_H2SC_X -1.880000e-01 // X of housing 2 [m]
#define EOM_H2SC_Y 0.000000e+00 // Y of housing 2 [m]
#define EOM_H2SC_Z 5.325000e-01 // Z of housing 2 [m]

// Spacecraft Mass Properties
#define EOM_SC_M 4.234930e+02 // SC mass [kg]
#define EOM_SC_IB_XX 1.891515e+01 // XX component of spacecraft moment of inertia tensor about B [kg m^(2)]
#define EOM_SC_IB_XY -1.368691e+00 // XY component of spacecraft moment of inertia tensor about B [kg m^(2)]
#define EOM_SC_IB_XZ -2.848753e+00 // XZ component of spacecraft moment of inertia tensor about B [kg m^(2)]
#define EOM_SC_IB_YX -1.368691e+00 // YX component of spacecraft moment of inertia tensor about B [kg m^(2)]
#define EOM_SC_IB_YY 1.876681e+01 // YY component of spacecraft moment of inertia tensor about B [kg m^(2)]
#define EOM_SC_IB_YZ -1.976778e+00 // YZ component of spacecraft moment of inertia tensor about B [kg m^(2)]
#define EOM_SC_IB_ZX -2.848753e+00 // ZX component of spacecraft moment of inertia tensor about B [kg m^(2)]
#define EOM_SC_IB_ZY -1.966778e+00 // ZY component of spacecraft moment of inertia tensor about B [kg m^(2)]
#define EOM_SC_IB_ZZ 2.011069e+02 // ZZ component of spacecraft moment of inertia tensor about B [kg m^(2)]
#define EOM_SC_IH1_XX 1.963625e+01 // XX component of spacecraft moment of inertia tensor about H1 [kg m^(2)]
#define EOM_SC_IH1_XY -1.181266e+00 // XY component of spacecraft moment of inertia tensor about H1 [kg m^(2)]
#define EOM_SC_IH1_XZ -6.105918e+00 // XZ component of spacecraft moment of inertia tensor about H1 [kg m^(2)]
#define EOM_SC_IH1_YX -1.181266e+00 // YX component of spacecraft moment of inertia tensor about H1 [kg m^(2)]
#define EOM_SC_IH1_YY 3.424660e+01 // YY component of spacecraft moment of inertia tensor about H1 [kg m^(2)]
#define EOM_SC_IH1_YZ -1.935420e+00 // YZ component of spacecraft moment of inertia tensor about H1 [kg m^(2)]
#define EOM_SC_IH1_ZX -6.105918e+00 // ZX component of spacecraft moment of inertia tensor about H1 [kg m^(2)]
#define EOM_SC_IH1_ZY -1.925420e+00 // ZY component of spacecraft moment of inertia tensor about H1 [kg m^(2)]
#define EOM_SC_IH1_ZZ 2.158703e+02 // ZZ component of spacecraft moment of inertia tensor about H1 [kg m^(2)]
#define EOM_SC_IH2_XX 1.963625e+01 // XX component of spacecraft moment of inertia tensor about H2 [kg m^(2)]
#define EOM_SC_IH2_XY -1.558734e+00 // XY component of spacecraft moment of inertia tensor about H2 [kg m^(2)]
#define EOM_SC_IH2_XZ 4.539004e-01 // XZ component of spacecraft moment of inertia tensor about H2 [kg m^(2)]
#define EOM_SC_IH2_YX -1.558734e+00 // YX component of spacecraft moment of inertia tensor about H2 [kg m^(2)]
#define EOM_SC_IH2_YY 3.466178e+01 // YY component of spacecraft moment of inertia tensor about H2 [kg m^(2)]
#define EOM_SC_IH2_YZ -1.935420e+00 // YZ component of spacecraft moment of inertia tensor about H2 [kg m^(2)]
#define EOM_SC_IH2_ZX 4.539004e-01 // ZX component of spacecraft moment of inertia tensor about H2 [kg m^(2)]
#define EOM_SC_IH2_ZY -1.925420e+00 // ZY component of spacecraft moment of inertia tensor about H2 [kg m^(2)]
#define EOM_SC_IH2_ZZ 2.162855e+02 // ZZ component of spacecraft moment of inertia tensor about H2 [kg m^(2)]

// Cold gas A thruster positions
#define COLDGAS_A_01_RX 8.109270e-01 // M-frame x coordinate for mount point of cold gas A thruster #01 [m]
#define COLDGAS_A_01_RY 6.392940e-01 // M-frame y coordinate for mount point of cold gas A thruster #01 [m]
#define COLDGAS_A_01_RZ 3.991360e-01 // M-frame z coordinate for mount point of cold gas A thruster #01 [m]
#define COLDGAS_A_02_RX 9.591110e-01 // M-frame x coordinate for mount point of cold gas A thruster #02 [m]
#define COLDGAS_A_02_RY 3.826360e-01 // M-frame y coordinate for mount point of cold gas A thruster #02 [m]
#define COLDGAS_A_02_RZ 3.991340e-01 // M-frame z coordinate for mount point of cold gas A thruster #02 [m]
#define COLDGAS_A_03_RX -9.591080e-01 // M-frame x coordinate for mount point of cold gas A thruster #03 [m]
#define COLDGAS_A_03_RY 3.826360e-01 // M-frame y coordinate for mount point of cold gas A thruster #03 [m]
#define COLDGAS_A_03_RZ 3.991360e-01 // M-frame z coordinate for mount point of cold gas A thruster #03 [m]
#define COLDGAS_A_04_RX -8.109280e-01 // M-frame x coordinate for mount point of cold gas A thruster #04 [m]
#define COLDGAS_A_04_RY 6.392970e-01 // M-frame y coordinate for mount point of cold gas A thruster #04 [m]
#define COLDGAS_A_04_RZ 3.991340e-01 // M-frame z coordinate for mount point of cold gas A thruster #04 [m]
#define COLDGAS_A_05_RX 1.481810e-01 // M-frame x coordinate for mount point of cold gas A thruster #05 [m]
#define COLDGAS_A_05_RY -1.021930e+00 // M-frame y coordinate for mount point of cold gas A thruster #05 [m]
#define COLDGAS_A_05_RZ 3.991360e-01 // M-frame z coordinate for mount point of cold gas A thruster #05 [m]
#define COLDGAS_A_06_RX -1.481830e-01 // M-frame x coordinate for mount point of cold gas A thruster #06 [m]
#define COLDGAS_A_06_RY -1.021932e+00 // M-frame y coordinate for mount point of cold gas A thruster #06 [m]
#define COLDGAS_A_06_RZ 3.991340e-01 // M-frame z coordinate for mount point of cold gas A thruster #06 [m]

// Cold gas A thruster directions
#define COLDGAS_A_01_ELV 5.198290e-01 // M-frame elevation (angle relative to xy-plane) for thrust vector of cold gas A thruster #01 [rad]
#define COLDGAS_A_01_AZM -1.661147e+00 // M-frame azimuth (angle relative to x-axis) for thrust vector of cold gas A thruster #01 [rad]
#define COLDGAS_A_02_ELV 5.198290e-01 // M-frame elevation (angle relative to xy-plane) for thrust vector of cold gas A thruster #02 [rad]
#define COLDGAS_A_02_AZM 2.708345e+00 // M-frame azimuth (angle relative to x-axis) for thrust vector of cold gas A thruster #02 [rad]
#define COLDGAS_A_03_ELV 5.198290e-01 // M-frame elevation (angle relative to xy-plane) for thrust vector of cold gas A thruster #03 [rad]
#define COLDGAS_A_03_AZM 4.332477e-01 // M-frame azimuth (angle relative to x-axis) for thrust vector of cold gas A thruster #03 [rad]
#define COLDGAS_A_04_ELV 5.198290e-01 // M-frame elevation (angle relative to xy-plane) for thrust vector of cold gas A thruster #04 [rad]
#define COLDGAS_A_04_AZM -1.480445e+00 // M-frame azimuth (angle relative to x-axis) for thrust vector of cold gas A thruster #04 [rad]
#define COLDGAS_A_05_ELV 5.198290e-01 // M-frame elevation (angle relative to xy-plane) for thrust vector of cold gas A thruster #05 [rad]
#define COLDGAS_A_05_AZM 2.527643e+00 // M-frame azimuth (angle relative to x-axis) for thrust vector of cold gas A thruster #05 [rad]
#define COLDGAS_A_06_ELV 5.198290e-01 // M-frame elevation (angle relative to xy-plane) for thrust vector of cold gas A thruster #06 [rad]
#define COLDGAS_A_06_AZM 6.139499e-01 // M-frame azimuth (angle relative to x-axis) for thrust vector of cold gas A thruster #06 [rad]

// Cold gas B thruster positions
#define COLDGAS_B_01_RX 8.125020e-01 // M-frame x coordinate for mount point of cold gas B thruster #01 [m]
#define COLDGAS_B_01_RY 6.365680e-01 // M-frame y coordinate for mount point of cold gas B thruster #01 [m]
#define COLDGAS_B_01_RZ 2.709170e-01 // M-frame z coordinate for mount point of cold gas B thruster #01 [m]
#define COLDGAS_B_02_RX 9.575370e-01 // M-frame x coordinate for mount point of cold gas B thruster #02 [m]
#define COLDGAS_B_02_RY 3.853620e-01 // M-frame y coordinate for mount point of cold gas B thruster #02 [m]
#define COLDGAS_B_02_RZ 2.709160e-01 // M-frame z coordinate for mount point of cold gas B thruster #02 [m]
#define COLDGAS_B_03_RX -9.575350e-01 // M-frame x coordinate for mount point of cold gas B thruster #03 [m]
#define COLDGAS_B_03_RY 3.853630e-01 // M-frame y coordinate for mount point of cold gas B thruster #03 [m]
#define COLDGAS_B_03_RZ 2.709170e-01 // M-frame z coordinate for mount point of cold gas B thruster #03 [m]
#define COLDGAS_B_04_RX -8.125020e-01 // M-frame x coordinate for mount point of cold gas B thruster #04 [m]
#define COLDGAS_B_04_RY 6.365700e-01 // M-frame y coordinate for mount point of cold gas B thruster #04 [m]
#define COLDGAS_B_04_RZ 2.709160e-01 // M-frame z coordinate for mount point of cold gas B thruster #04 [m]
#define COLDGAS_B_05_RX 1.450340e-01 // M-frame x coordinate for mount point of cold gas B thruster #05 [m]
#define COLDGAS_B_05_RY -1.021930e+00 // M-frame y coordinate for mount point of cold gas B thruster #05 [m]
#define COLDGAS_B_05_RZ 2.709170e-01 // M-frame z coordinate for mount point of cold gas B thruster #05 [m]
#define COLDGAS_B_06_RX -1.450350e-01 // M-frame x coordinate for mount point of cold gas B thruster #06 [m]
#define COLDGAS_B_06_RY -1.021932e+00 // M-frame y coordinate for mount point of cold gas B thruster #06 [m]
#define COLDGAS_B_06_RZ 2.709160e-01 // M-frame z coordinate for mount point of cold gas B thruster #06 [m]

// Cold gas B thruster directions
#define COLDGAS_B_01_ELV 5.198293e-01 // M-frame elevation (angle relative to xy-plane) for thrust vector of cold gas B thruster #01 [rad]
#define COLDGAS_B_01_AZM -1.661147e+00 // M-frame azimuth (angle relative to x-axis) for thrust vector of cold gas B thruster #01 [rad]
#define COLDGAS_B_02_ELV 5.198293e-01 // M-frame elevation (angle relative to xy-plane) for thrust vector of cold gas B thruster #02 [rad]
#define COLDGAS_B_02_AZM 2.708345e+00 // M-frame azimuth (angle relative to x-axis) for thrust vector of cold gas B thruster #02 [rad]
#define COLDGAS_B_03_ELV 5.198293e-01 // M-frame elevation (angle relative to xy-plane) for thrust vector of cold gas B thruster #03 [rad]
#define COLDGAS_B_03_AZM 4.332476e-01 // M-frame azimuth (angle relative to x-axis) for thrust vector of cold gas B thruster #03 [rad]
#define COLDGAS_B_04_ELV 5.198293e-01 // M-frame elevation (angle relative to xy-plane) for thrust vector of cold gas B thruster #04 [rad]
#define COLDGAS_B_04_AZM -1.480445e+00 // M-frame azimuth (angle relative to x-axis) for thrust vector of cold gas B thruster #04 [rad]
#define COLDGAS_B_05_ELV 5.198293e-01 // M-frame elevation (angle relative to xy-plane) for thrust vector of cold gas B thruster #05 [rad]
#define COLDGAS_B_05_AZM 2.527642e+00 // M-frame azimuth (angle relative to x-axis) for thrust vector of cold gas B thruster #05 [rad]
#define COLDGAS_B_06_ELV 5.198293e-01 // M-frame elevation (angle relative to xy-plane) for thrust vector of cold gas B thruster #06 [rad]
#define COLDGAS_B_06_AZM 6.139502e-01 // M-frame azimuth (angle relative to x-axis) for thrust vector of cold gas B thruster #06 [rad]

// CMNT thruster positions
#define CMNT_01_RX 1.012104e+00 // M-frame x coordinate for mount point of cold gas thruster #01 [m]
#define CMNT_01_RY 1.150400e-01 // M-frame y coordinate for mount point of cold gas thruster #01 [m]
#define CMNT_01_RZ 6.027550e-01 // M-frame z coordinate for mount point of cold gas thruster #01 [m]
#define CMNT_02_RX 1.016708e+00 // M-frame x coordinate for mount point of cold gas thruster #02 [m]
#define CMNT_02_RY 8.798800e-02 // M-frame y coordinate for mount point of cold gas thruster #02 [m]
#define CMNT_02_RZ 3.197190e-01 // M-frame z coordinate for mount point of cold gas thruster #02 [m]
#define CMNT_03_RX 1.016752e+00 // M-frame x coordinate for mount point of cold gas thruster #03 [m]
#define CMNT_03_RY -8.798800e-02 // M-frame y coordinate for mount point of cold gas thruster #03 [m]
#define CMNT_03_RZ 3.197190e-01 // M-frame z coordinate for mount point of cold gas thruster #03 [m]
#define CMNT_04_RX 1.013120e+00 // M-frame x coordinate for mount point of cold gas thruster #04 [m]
#define CMNT_04_RY 8.798800e-02 // M-frame y coordinate for mount point of cold gas thruster #04 [m]
#define CMNT_04_RZ 6.027550e-01 // M-frame z coordinate for mount point of cold gas thruster #04 [m]
#define CMNT_05_RX -1.012104e+00 // M-frame x coordinate for mount point of cold gas thruster #05 [m]
#define CMNT_05_RY -1.150400e-01 // M-frame y coordinate for mount point of cold gas thruster #05 [m]
#define CMNT_05_RZ 6.027550e-01 // M-frame z coordinate for mount point of cold gas thruster #05 [m]
#define CMNT_06_RX -1.016708e+00 // M-frame x coordinate for mount point of cold gas thruster #06 [m]
#define CMNT_06_RY -8.798800e-02 // M-frame y coordinate for mount point of cold gas thruster #06 [m]
#define CMNT_06_RZ 3.197190e-01 // M-frame z coordinate for mount point of cold gas thruster #06 [m]
#define CMNT_07_RX -1.016752e+00 // M-frame x coordinate for mount point of cold gas thruster #07 [m]
#define CMNT_07_RY 8.798800e-02 // M-frame y coordinate for mount point of cold gas thruster #07 [m]
#define CMNT_07_RZ 3.197190e-01 // M-frame z coordinate for mount point of cold gas thruster #07 [m]
#define CMNT_08_RX -1.086752e+00 // M-frame x coordinate for mount point of cold gas thruster #08 [m]
#define CMNT_08_RY 1.150400e-01 // M-frame y coordinate for mount point of cold gas thruster #08 [m]
#define CMNT_08_RZ 6.027550e-01 // M-frame z coordinate for mount point of cold gas thruster #08 [m]

// CMNT thruster directions
#define CMNT_01_ELV 3.608280e-01 // M-frame elevation (angle relative to xy-plane) for thrust vector of cold gas thruster #01 [rad]
#define CMNT_01_AZM 9.428325e-01 // M-frame azimuth (angle relative to x-axis) for thrust vector of cmnt thruster #01 [rad]
#define CMNT_02_ELV -6.177813e-01 // M-frame elevation (angle relative to xy-plane) for thrust vector of cold gas thruster #02 [rad]
#define CMNT_02_AZM 7.903010e-01 // M-frame azimuth (angle relative to x-axis) for thrust vector of cmnt thruster #02 [rad]
#define CMNT_03_ELV -6.177813e-01 // M-frame elevation (angle relative to xy-plane) for thrust vector of cold gas thruster #03 [rad]
#define CMNT_03_AZM -7.904010e-01 // M-frame azimuth (angle relative to x-axis) for thrust vector of cmnt thruster #03 [rad]
#define CMNT_04_ELV 3.608280e-01 // M-frame elevation (angle relative to xy-plane) for thrust vector of cold gas thruster #04 [rad]
#define CMNT_04_AZM -9.428325e-01 // M-frame azimuth (angle relative to x-axis) for thrust vector of cmnt thruster #04 [rad]
#define CMNT_05_ELV 3.608280e-01 // M-frame elevation (angle relative to xy-plane) for thrust vector of cold gas thruster #05 [rad]
#define CMNT_05_AZM 2.198760e+00 // M-frame azimuth (angle relative to x-axis) for thrust vector of cmnt thruster #05 [rad]
#define CMNT_06_ELV -6.177813e-01 // M-frame elevation (angle relative to xy-plane) for thrust vector of cold gas thruster #06 [rad]
#define CMNT_06_AZM -2.351292e+00 // M-frame azimuth (angle relative to x-axis) for thrust vector of cmnt thruster #06 [rad]
#define CMNT_07_ELV -6.177813e-01 // M-frame elevation (angle relative to xy-plane) for thrust vector of cold gas thruster #07 [rad]
#define CMNT_07_AZM 2.351292e+00 // M-frame azimuth (angle relative to x-axis) for thrust vector of cmnt thruster #07 [rad]
#define CMNT_08_ELV -6.177813e-01 // M-frame elevation (angle relative to xy-plane) for thrust vector of cold gas thruster #08 [rad]
#define CMNT_08_AZM 2.351292e+00 // M-frame azimuth (angle relative to x-axis) for thrust vector of cmnt thruster #08 [rad]

// Typical Noise levels
#define NOISE_COLD_GAS 1.000000e-07 // Amplitude spectral density of force noise for a typical single Cold Gas Thruster [kg m s^(-2) Hz^(-1/2)]
#define NOISE_CMNT 1.000000e-07 // Amplitude spectral density of force noise for a typical single CMNT Thruster [kg m s^(-2) Hz^(-1/2)]
#define NOISE_GRS_POS 3.500000e-09 // Amplitude spectral density of displacement noise for GRS capacitive sensing in linear DoFs [m Hz^(-1/2)]
#define NOISE_GRS_ANG 2.500000e-07 // Amplitude spectral density of displacement noise for GRS capacitive sensing in angular DoFs [rad Hz^(-1/2)]
#define NOISE_OMS_POS 1.000000e-12 // Amplitude spectral density of displacement noise for OMS optical sensing in linear DoFs [m Hz^(-1/2)]
#define NOISE_OMS_ANG 1.000000e-08 // Amplitude spectral density of displacement noise for OMS optical sensing in angular DoFs [rad Hz^(-1/2)]
