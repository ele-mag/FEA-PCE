import numpy as np
import math
import pandas as pd
sin = math.sin
cos = math.cos
log = math.log
pi = math.pi
import os


def calculate_jianhe(gonggao,kuaju,zhijing,genshu,jianju,f_start, fstop, fpoints,file_path):
    """
    model
    """
    # f = 12e9
    c0 = 3e8  # c
    u0 = 4e-7 * pi
    e0 = 8.854e-12

    w1 = 244e-6
    w2 = 244e-6
    s = 0.29 * w1
    s=jianju/10**3

    p1 = 56e-6
    p2 = 56e-6
    p1=p2=kuaju/8/10**3
    # print(p1)
    g = 224e-6
    g=kuaju / 5*2 / 10 ** 3
   # print(g)
    d1 = 0e-6
    d2 = 0e-6

    d_b = p1 + p2 + d1 + d2 + g
    h1 = 254e-6
    h2 = 300e-6
    hb = 344e-6
    hb =gonggao/10**3
    rw = 8.5e-6
    rw = zhijing/10**3/2
    er1 = 9.86
    er2 = 9.86
    Z0 = 50

    Freq = []
    S11_list = []
    S21_list = []
    f = np.arange(5e9, 40e9 + 0.1e9, 0.1e9)

    directory = os.path.dirname(file_path)


    if not os.path.exists(directory):
        os.makedirs(directory)

    with open(file_path, "w") as file:
        file.write("# GHz S RI R 50\n")


        for f in f:

            beita_0 = 2 * pi * f * (u0 * e0) ** 0.5

            a = 1 - (hb - h2) / (hb - h1)
            b = -2 * d_b
            c = d_b * d_b + (h1 - h2) * (hb - h2)
            if a == 0:
                x_c = d_b / 2
            else:  # b*b-4*a*c >= 0
                x_c = -(b + (b * b - 4 * a * c) ** 0.5) / (2 * a)

            y_c = (hb * hb - h1 * h1 - x_c * x_c) / (2 * (hb - h1))
            r_c = hb - y_c

            hx_0 = y_c + (r_c * r_c - (0 - x_c) * (0 - x_c)) ** 0.5

            x_BC = p1 + d1 / 2
            x_CD = (x_c + p1 + d1) / 2
            x_DE = (x_c + p1 + d1 + g) / 2
            x_EF = p1 + d1 + g + d2 / 2

            # hx = y_c+(r_c*r_c-(x-x_c)*(x-x_c))**0.5

            hx_BC = y_c + (r_c * r_c - (x_BC - x_c) * (x_BC - x_c)) ** 0.5
            hx_CD = y_c + (r_c * r_c - (x_CD - x_c) * (x_CD - x_c)) ** 0.5
            hx_DE = y_c + (r_c * r_c - (x_DE - x_c) * (x_DE - x_c)) ** 0.5
            hx_EF = y_c + (r_c * r_c - (x_EF - x_c) * (x_EF - x_c)) ** 0.5
            h_w1 = y_c + (r_c * r_c - (p1 / 2 - x_c) * (p1 / 2 - x_c)) ** 0.5 - h1
            h_w2 = y_c + (r_c * r_c - ((d_b - p2 / 2) - x_c) * ((d_b - p2 / 2) - x_c)) ** 0.5

            # def z_calculation(h):
            #     u0 = log(h/rw+((h/rw)**2-1)**0.5)
            #     k0 = log((((h/rw)**2-1)**0.5+hs/rw)/(((h/rw)**2-1)**0.5-hs/rw))
            #     ew = (1-(k0/u0)*((er-1)/er))

            er_1 = (er1 * h1 + 1 * (hx_BC - h1)) / hx_BC

            u0_BC = log(hx_BC / rw + ((hx_BC / rw) ** 2 - 1) ** 0.5)
            k0_BC = log((((hx_BC / rw) ** 2 - 1) ** 0.5 + h1 / rw) / (((hx_BC / rw) ** 2 - 1) ** 0.5 - h1 / rw))
            ew_BC = (1 - (k0_BC / u0_BC) * ((er_1 - 1) / er_1)) ** (-0.5)
            Zc_BC = 30 * log((2 * hx_BC * (s * s + 4 * hx_BC * hx_BC) ** 0.5) / (s * rw))

            theta_BC = beita_0 * d1 * ew_BC ** 0.5
            A_BC = np.array(
                [[cos(theta_BC), complex(0, Zc_BC * sin(theta_BC))], [complex(0, (sin(theta_BC)) / Zc_BC), cos(theta_BC)]])

            u0_CD = log(hx_CD / rw + ((hx_CD / rw) ** 2 - 1) ** 0.5)
            k0_CD = log((((hx_CD / rw) ** 2 - 1) ** 0.5 + 0 / rw) / (((hx_CD / rw) ** 2 - 1) ** 0.5 - 0 / rw))
            ew_CD = (1 - (k0_CD / u0_CD) * ((1 - 1) / 1))
            Zc_CD = 30 * log((2 * hx_CD * (s * s + 4 * hx_CD * hx_CD) ** 0.5) / (s * rw))

            theta_CD = beita_0 * (x_c - p1 - d1)

            A_CD = np.array(
                [[cos(theta_CD), complex(0, Zc_CD * sin(theta_CD))], [complex(0, (sin(theta_CD)) / Zc_CD), cos(theta_CD)]])

            u0_DE = log(hx_DE / rw + ((hx_DE / rw) ** 2 - 1) ** 0.5)
            k0_DE = log((((hx_DE / rw) ** 2 - 1) ** 0.5 + 0 / rw) / (((hx_DE / rw) ** 2 - 1) ** 0.5 - 0 / rw))
            ew_DE = (1 - (k0_DE / u0_DE) * ((1 - 1) / 1))
            Zc_DE = 30 * log((2 * hx_DE * (s * s + 4 * hx_DE * hx_DE) ** 0.5) / (s * rw))

            theta_DE = beita_0 * (p1 + d1 + g - x_c)

            A_DE = np.array(
                [[cos(theta_DE), complex(0, Zc_DE * sin(theta_DE))], [complex(0, (sin(theta_DE)) / Zc_DE), cos(theta_DE)]])

            er_2 = (er2 * h2 + 1 * (hx_EF - h2)) / hx_EF

            u0_EF = log(hx_EF / rw + ((hx_EF / rw) ** 2 - 1) ** 0.5)
            k0_EF = log((((hx_EF / rw) ** 2 - 1) ** 0.5 + h2 / rw) / (((hx_EF / rw) ** 2 - 1) ** 0.5 - h2 / rw))
            ew_EF = (1 - (k0_EF / u0_EF) * ((er_2 - 1) / er_2))
            Zc_EF = 30 * log((2 * hx_EF * (s * s + 4 * hx_EF * hx_EF) ** 0.5) / (s * rw))

            theta_EF = beita_0 * d2 * ew_EF ** 0.5

            A_EF = np.array(
                [[cos(theta_EF), complex(0, Zc_EF * sin(theta_EF))], [complex(0, (sin(theta_EF)) / Zc_EF), cos(theta_EF)]])


            # e_1 = 6.636  #
            # L_s1 = 0.5*e_1**0.5*p1/c0*Z0
            # L_w1 = L_s1 + u0/(2*pi)*p1*log((y_c+(r_c*r_c-(p1/2-x_c)*(p1/2-x_c))**0.5-h1)/rw+((y_c+(r_c*r_c-(p1/2-x_c)*(p1/2-x_c))**0.5-h1)**2-1)**0.5)
            # C_1 = e_1**0.5*p1/(c0*Z0)

            e_1 = 6.636
            Ls1 = 0.5 * e_1 ** 0.5 * p1 * Z0 / c0
            Lw1 = 0.5 * e_1 ** 0.5 * p1 * Z0 / c0
            C1 = e_1 ** 0.5 * p1 / (c0 * Z0) + 7.5e-15

            Z_Ls1 = complex(0, 2 * pi * f * Ls1)
            Z_Lw1 = complex(0, 2 * pi * f * Lw1) + (u0 * p1 / (4 * pi)) * log(
                (2 * h_w1 * (s * s + 4 * h_w1 * h_w1) ** 0.5) / (s * rw))
            Y_C1 = complex(0, 2 * pi * f * C1)

            A_AB = np.array([[1, Z_Ls1], [0, 1]]) @ np.array([[1, 0], [Y_C1, 1]]) @ np.array([[1, Z_Lw1], [0, 1]])

            e_2 = 6.636
            Ls2 = 0.5 * e_2 ** 0.5 * p2 * Z0 / c0
            Lw2 = 0.5 * e_2 ** 0.5 * p2 * Z0 / c0 + (u0 * p1 / (4 * pi)) * log(
                (2 * h_w2 * (s * s + 4 * h_w2 * h_w2) ** 0.5) / (s * rw))
            C2 = e_2 ** 0.5 * p2 / (c0 * Z0) + 7.5e-15

            Z_Ls2 = complex(0, 2 * pi * f * Ls2)
            Z_Lw2 = complex(0, 2 * pi * f * Lw2)
            Y_C2 = complex(0, 2 * pi * f * C2)

            A_FG = np.array([[1, Z_Lw2], [0, 1]]) @ np.array([[1, 0], [Y_C2, 1]]) @ np.array([[1, Z_Ls2], [0, 1]])

            A_Total = A_AB @ A_BC @ A_CD @ A_DE @ A_EF @ A_FG
            # print(A_Total)
            AA = A_Total[0, 0]
            BB = A_Total[0, 1]
            CC = A_Total[1, 0]
            DD = A_Total[1, 1]
            S11_Complex = (AA * Z0 + BB - CC * Z0 * Z0 - DD * Z0) / (AA * Z0 + BB + CC * Z0 * Z0 + DD * Z0)
            S21_Complex = (2 * Z0) / (AA * Z0 + BB + CC * Z0 * Z0 + DD * Z0)
            S12_Complex = (2 * (AA * DD - BB * CC) * Z0) / (AA * Z0 + BB + CC * Z0 * Z0 + DD * Z0)
            S22_Complex = (-AA * Z0 + BB - CC * Z0 * Z0 + DD * Z0) / (AA * Z0 + BB + CC * Z0 * Z0 + DD * Z0)
            f_in_ghz = f / 1e9
            file.write(f"{f_in_ghz:.9f} "
                       f"{S11_Complex.real:.9e} {S11_Complex.imag:.9e} "
                       f"{S12_Complex.real:.9e} {S12_Complex.imag:.9e} "
                       f"{S21_Complex.real:.9e} {S21_Complex.imag:.9e} "
                       f"{S22_Complex.real:.9e} {S22_Complex.imag:.9e} "
                       f"\n")


            # print("freq=",f/1e9,"Ghz")
            # print("S11=", S11_dB, "dB")
            # print("S21=", S21_dB, "dB")
            # print("S12=", S12_dB, "dB")
            # print("S22=", S22_dB, "dB")
            # S21_Phase = np.angle(S21_Complex,deg=True)
            # print(S21_Phase)
        # print(Freq)
        # print(S11_list)
        # print(S21_list)

    # df = pd.DataFrame({'Freq/GHz':Freq,'S11/dB':S11_list,'S21/dB':S21_list})
    # df.to_excel('S_Parameter.xlsx',index=False)

# calculate_jianhe(0.25, 0.6, 0.025, 2, 0.11,freq_data["start"], freq_data["stop"], freq_data["points"])

