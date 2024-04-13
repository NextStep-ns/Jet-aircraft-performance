#FLIGHT MECHANICS

from math import sqrt, asin, pi, log, cos, radians
import colorama
from colorama import Fore
import matplotlib.pyplot as plt
from Formulas import *
import numpy
import pandas as pd
# import matplotlib.animation as animation

################################# TAKE-OFF #####################"
print(Fore.RED + '                          TAKE-OFF                       '+Fore.RESET)
print('')

def Rho(altitude):
    return Pressure(altitude)/(287.05*Temperature(altitude))
print('rho=', Rho(Feet_Into_Meter(1900)))

def Pitch(alpha, delta_e):
    return 0.24 - 0.18*alpha + 0.28*delta_e

def Thrust(delta_e, rho):
    return delta_e*((rho/1.225)**0.6)*12980
print('thrust=', Thrust(1, Rho(Feet_Into_Meter(1900))))

def V_stall(lift_coefficient, rho, mass):
    return sqrt((2*mass*g)/(rho*Wing_area*lift_coefficient))
print('V_stall=',V_stall(Max_Lift_Coeff_Flap_Partially_Extended_TO, Rho(Feet_Into_Meter(1900)), Max_weight))

def V2(lift_coefficient, rho, mass):
    return 1.2*V_stall(lift_coefficient, rho, mass)
print('V2=', V2(Max_Lift_Coeff_Flap_Partially_Extended_TO, Rho(Feet_Into_Meter(1900)), Max_weight))

def Cl_TO(rho, velocity):
    return (2*Max_weight*g)/((rho)*Wing_area*velocity**2)
print('CL=', Cl_TO(Rho(Feet_Into_Meter(1900)), V2(Max_Lift_Coeff_Flap_Partially_Extended_TO, Rho(Feet_Into_Meter(1900)), Max_weight)))

def Cd_TO(Cl):
    return Polar(Cl)
print('CD=', Cd_TO(Cl_TO(Rho(Feet_Into_Meter(1900)), V2(Max_Lift_Coeff_Flap_Partially_Extended_TO, Rho(Feet_Into_Meter(1900)), Max_weight))))

def CL_CD(rho, velocity, Cl):
    return Cl_TO(rho, velocity)/Cd_TO(Cl)
print('CL/CD = ', CL_CD(Rho(Feet_Into_Meter(1900)), V2(Max_Lift_Coeff_Flap_Partially_Extended_TO, Rho(Feet_Into_Meter(1900)), Max_weight), Cl_TO(Rho(Feet_Into_Meter(1900)), V2(Max_Lift_Coeff_Flap_Partially_Extended_TO, Rho(Feet_Into_Meter(1900)), Max_weight))))

############################## RESULTS ##############################

def Lenght_Of_Flight():
    return (10.5+(((V2(Max_Lift_Coeff_Flap_Partially_Extended_TO, Rho(Feet_Into_Meter(1900)), Max_weight)**2)-(V_lof**2)))/(2*g))/(((Thrust(1, Rho(Feet_Into_Meter(1900)))/((Max_weight*g)))-(1/CL_CD(Rho(Feet_Into_Meter(1900)), V2(Max_Lift_Coeff_Flap_Partially_Extended_TO, Rho(Feet_Into_Meter(1900)), Max_weight), Cl_TO(Rho(Feet_Into_Meter(1900)), V2(Max_Lift_Coeff_Flap_Partially_Extended_TO, Rho(Feet_Into_Meter(1900)), Max_weight))))))

def Length_Of_Rolling():
    return (V_lof**2)/(g*(((Thrust(1, Rho(Feet_Into_Meter(1900))))/(Max_weight*g))-r))

def Take_Off_Distance():
    return Length_Of_Rolling()+Lenght_Of_Flight()

distance_dict = {'Length of flight': Lenght_Of_Flight(),
                 'Length of rolling': Length_Of_Rolling(),
                 'Take off distance': Take_Off_Distance()}

distance = pd.Series(distance_dict)

blabla= pd.DataFrame({'distance': distance})
print(blabla)
distance.to_csv(r'C:/Users/Lenovo/Documents/works/3RD YEAR/1st semester/Fligh Mechanics/DM/CSV.csv')


print(Fore.RED +'                          CLIMBING                             ' + Fore.WHITE)
print('    ')
print(Fore.YELLOW +'Question #1'+ Fore.RESET)
print('    ')

######################################################### CLIMBING- Q1 #################################################

def Rho_FL350():
    return Rho(Feet_Into_Meter(35000))

Rho_average =(Rho_FL350()+Rho(Feet_Into_Meter(1900)))/2

def Cl2(rho, velocity):
    return (2*Max_weight*g)/((rho)*Wing_area*(velocity**2))

def Cd2():
    return Polar(Cl2(Rho_average, 88))

def Sin_Gamma():
    return vario/v_climbing

def Thrust_required(rho, sin_gamma, Cd, velocity):
    return 0.5*(rho)*Wing_area*(velocity**2)*Cd + Max_weight*g*sin_gamma
print(Thrust_required(Rho_average, Sin_Gamma(), Cd2(), 88))

rho_average_dict = {'Rho_average': Rho_average}
rho_average = pd.Series(rho_average_dict)
Cl2_dict = {'Cl2': Cl2(Rho_average, 88),}
cl2 = pd.Series(Cl2_dict)
Cd2_dict = {'Cd2': Cd2()}
cd2 = pd.Series(Cd2_dict)
sin_gamma_dict = {'sin gamma': Sin_Gamma()}
sin_gamma = pd.Series(sin_gamma_dict)
Thrust_required_dict = {'Thrust required': Thrust_required(Rho_average, Sin_Gamma(), Cd2(), 88)}
thrust_required = pd.Series(Thrust_required_dict)


question1_CLIMBING= pd.DataFrame({'rho': rho_average, 'Cl': cl2, 'Cd': cd2, 'sinus gamma': sin_gamma, 'Thrust': thrust_required})
print(question1_CLIMBING)

question1_CLIMBING.to_csv(r'C:/Users/Lenovo/Documents/works/3RD YEAR/1st semester/Fligh Mechanics/DM/CSV.csv',sep =  ',')

############################# CLIMBING-Q2 ####################################

print('  ')
print(Fore.YELLOW +'Question #2'+ Fore.RESET)
print('    ')

Thrust_max = Thrust(1, (Rho_FL350()+Rho(Feet_Into_Meter(1900)))/2)

def Sin_Gamma_For_Max_Thrust(rho, Cd, thrust, velocity):
    return (-0.5*rho*Wing_area*(velocity**2)*Cd+thrust)/(Max_weight*g)

Resulting_slope = asin(Sin_Gamma_For_Max_Thrust(Rho_average, Cd2(), Thrust_max, 88))*(180/pi)

def Rate_Of_Climb():
    return Meter_Into_Feet(Sin_Gamma_For_Max_Thrust(Rho_average, Cd2(), Thrust_max, 88)*88)*60

max_thrust_dict = {'max_thrust': Thrust_max}
max_thrust = pd.Series(max_thrust_dict)
slope_dict = {'resulting_slope': Resulting_slope}
slope = pd.Series(slope_dict)
rate_of_climb_dict = {'rate of climb': Rate_Of_Climb()}
rate_of_climb = pd.Series(rate_of_climb_dict)

question2_CLIMBING= pd.DataFrame({'thrust': max_thrust, 'Slope': slope, 'Rate of climb': rate_of_climb})
print(question2_CLIMBING)

############################ CLIMBING-Q3 #####################################
print('')
print(Fore.YELLOW +'Question #3'+ Fore.RESET)
print('')

Rho_average_after_3000ft = (Rho(Feet_Into_Meter(3000))+Rho_FL350())/2

Cl3 = Cl2(Rho_average_after_3000ft, 65)

Cd3 = Polar(Cl3)
print(Cd3)

Thrust_max_with_1_engine = Thrust(1, Rho_average_after_3000ft)/2
print(Thrust_max_with_1_engine)

new_sin_gamma = Sin_Gamma_For_Max_Thrust(Rho_average_after_3000ft, Cd3, Thrust_max_with_1_engine, 65)

Resulting_slope2 = asin(new_sin_gamma)*(180/pi)
print(Resulting_slope2)

rho_average_3000ft_dict = {'rho average 3000ft': Rho_average_after_3000ft}
rho_average_3000ft = pd.Series(rho_average_3000ft_dict)
Cl3_dict = {'Cl3': Cl3}
cl3 = pd.Series(Cl3_dict)
Cd3_dict = {'Cd3': Cd3}
cd3 = pd.Series(Cd3_dict)
thrust_max_with_1_engine_dict = {'thrust_max_with_1_engine': Thrust_max_with_1_engine}
thrust_max_with_1_engine = pd.Series(thrust_max_with_1_engine_dict)
sin_gamma_1_engine_dict = {'sin_gamma_1_engine': new_sin_gamma}
sin_gamma_1_engine = pd.Series(sin_gamma_1_engine_dict)
slope_1_engine_dict = {'slope_1_engine': Resulting_slope2}
slope_1_engine = pd.Series(slope_1_engine_dict)

question3_CLIMBING= pd.DataFrame({'rho average': round(rho_average_3000ft,3),'thrust': round(thrust_max_with_1_engine, 3),'Lift coeff3': round(cl3, 3)  ,'Drag coeff': round(cd3, 3),'Slope with 1 engine': round(slope_1_engine, 3), 'Sin gamma 1 engine': round(sin_gamma_1_engine, 3)})
print(question3_CLIMBING)
##################################################### CRUISING FLIGHT ##################################################

########################### CRUISING FLIGHT - Q1 ########################
print('    ')
print(Fore.RED + '                    CRUISING FLIGHT                   ' + Fore.WHITE)
print('    ')
print(Fore.YELLOW +'Question #1'+ Fore.RESET)

def CLVmax(Fmax, mass):
    return ((1/(2*epsilon))*((Fmax/(mass*g))-sqrt(((Fmax/(mass*g))**2)-(4*epsilon*Cd0))))

def Vmax(mass, rho, CL):
    return sqrt((2*mass*g)/(rho*Wing_area*CL))

def Lifting_Ceiling():
    Ps = (2*Max_weight*g)/(gamma_lifting_ceiling*(mach_number**2)*Wing_area*Max_Lift_Coeff_Flap_Retracted)
    return (1-((Ps/101325))**(1/5.226))/(22.557*(10**-6))
Lifting_Ceiling()
print('Lifting_Ceiling : ',Lifting_Ceiling())

lifting_ceiling = numpy.linspace(Lifting_Ceiling(),Lifting_Ceiling(), 27)

Vstall_speed = []
altitude = []
def V_Stall_At_Each_Altitude():

    for i in range(0, int(Lifting_Ceiling()), 500):
        altitude.append(i)
        Vstall_speed.append(V_stall(Max_Lift_Coeff_Flap_Retracted, Rho(i), Max_weight))
    return Vstall_speed, altitude
V_Stall_At_Each_Altitude()


Max_speed = []
def V_Max_At_Each_Altitude():

    for i in range(0,int(Lifting_Ceiling()), 500):
        Fmax = Thrust(1, Rho(i))
        clvmax = CLVmax(Fmax, Max_weight)
        vmax = Vmax(Max_weight, Rho(i), clvmax)
        Max_speed.append(vmax)
    print('length of the list max speed :', len(Max_speed))
    return Max_speed
V_Max_At_Each_Altitude()

velocity = numpy.linspace(Vstall_speed[1], Max_speed[-1], 27)

########### GRAPH ###########

plt.grid(True)
plt.plot(Max_speed,altitude, "r", linewidth=0.8, marker="*", label="Trajet 1")
plt.plot(Vstall_speed,altitude, "b", linewidth=0.8, marker="*", label="Trajet 2")
plt.plot(velocity,lifting_ceiling, "g", linewidth=0.8, label="Trajet 3")
plt.xlabel('Velocity in m/s')
plt.ylabel('Altitude in meter')
plt.show()

################# STEADY CRUISE - Q3 ############################

print('    ')
print(Fore.YELLOW +'Question #3'+ Fore.RESET)

thrust_rq = []
def Thrust_Required_At_Each_Altitude():

    for i in range(0, int(Lifting_Ceiling()), 500):
        cl = Cl_TO(Rho(i), 165)
        cd = Cd_TO(cl)
        thrust_rq.append(Thrust_required(Rho(i), 0, cd, 165))
    return thrust_rq
Thrust_Required_At_Each_Altitude()

###### GRAPH ######
plt.grid(True)
plt.plot(altitude,thrust_rq, "g", linewidth=0.8, label="Trajet 1")
plt.xlabel('Altitude in meter')
plt.ylabel('Thrust in newton')
plt.show()
print(Fore.WHITE +'GRAPH'+ Fore.RESET)


################# STEADY CRUISE - Q4 ############################
print('    ')
print(Fore.YELLOW +'Question #4'+ Fore.RESET)


def Stamina(cs, cl, cd, mass):
    return (1/g)*(1/cs)*(cl/cd)*log((mass*g)/((mass-Fuel_Capacity)*g))

m_average = ((Max_weight + (Max_weight-Fuel_Capacity)))/2
Fmax_at_10000 = Thrust(1, Rho(10000))
clvmax_at_10000 = CLVmax(Fmax_at_10000, m_average)
vmax_at_10000 = Vmax(m_average,Rho(10000), clvmax_at_10000)

velocity = []
stamina_list = []
def Stamina_Graph():

    for i in range(1, int(vmax_at_10000)):
        velocity.append(i)
        cl_at_10000 = Cl_TO(Rho(10000), i)
        cd_at_10000 = Cd_TO(cl_at_10000)
        stamina_list.append(Stamina(cs, cl_at_10000, cd_at_10000, m_average)/3600)
    return stamina_list, velocity
Stamina_Graph()

######## GRAPH ########

plt.grid(True)
plt.plot(velocity,stamina_list, "g", linewidth=0.8, label="Trajet 1")
plt.xlabel('Velocity in m/s')
plt.ylabel('stamina in hour')
plt.show()

print('GRAPH')

################# STEADY CRUISE - Q5 #############################
print('    ')
print(Fore.YELLOW +'Question #5'+ Fore.RESET)

def Range(rho, cs, cl, cd, mass):
    return 2*(1/g)*sqrt((2/rho*Wing_area))*(1/cs)*(sqrt(cl)/cd)*(sqrt(mass*g)-sqrt((mass-Fuel_Capacity)*g))

range_list =[]
def Range_Graph():

    for i in range(1, int(vmax_at_10000)):
        range_list.append(Range(Rho(10000), cs, Cl_TO(Rho(10000), i), Cd_TO(Cl_TO(Rho(10000), i)), m_average))
    return range_list
Range_Graph()

plt.grid(True)
plt.plot(velocity,range_list, "g", linewidth=0.8, label="Trajet 1")
plt.xlabel('Velocity in m/s')
plt.ylabel('range')
plt.show()
print('GRAPH')


################################################# TURNING ##############################################################

print('    ')
print(Fore.RED +'                          TURNING                         '+ Fore.RESET)
print(Fore.YELLOW +'QUESTION 1'+ Fore.RESET)


def n(gamma):
    return (1/(cos(radians(gamma))))

def Vequivalent(n, velocity):
    return sqrt(abs(n))*velocity

def Thrust_turning(rho, veq, cd):
    return 0.5*rho*Wing_area*(veq**2)*cd

angle_list = []
required_thrust_with_inclination = []
def Required_Thrust_With_Inclination():

    cd_at_4500 = Cd_TO(Cl_TO(Rho(4500), 81))

    for i in range(1, 60):
        angle_list.append(i)
        required_thrust_with_inclination.append(Thrust_turning(Rho(4500), Vequivalent(n(i), 81), cd_at_4500))
    return required_thrust_with_inclination, angle_list

Required_Thrust_With_Inclination()

plt.grid(True)
plt.plot(angle_list,required_thrust_with_inclination, "g", linewidth=0.8, label="Trajet 1")
plt.xlabel('Angle in degree')
plt.ylabel('Thrust in newton')
plt.show()
print('GRAPH')

####################### TURNING - Q2 #############################

Vstall_LANDING = []
Vstall_TO = []
Vstall_rectracted = []
angle2 = []
def Stall_Speed_With_Flap_CO():
    for i in range(1, 60):
        angle2.append(i)
        Vstall_rectracted.append(n(i)*V_stall(Max_Lift_Coeff_Flap_Retracted, Rho(4500), Max_weight))
        Vstall_TO.append(n(i) * V_stall(Max_Lift_Coeff_Flap_Partially_Extended_TO, Rho(4500), Max_weight))
        Vstall_LANDING.append(n(i) * V_stall(Max_Lift_Coeff_Flap_Partially_Extended_LANDING, Rho(4500), Max_weight))
    return Vstall_rectracted,Vstall_TO,Vstall_LANDING, angle2
Stall_Speed_With_Flap_CO()

plt.grid(True)
plt.plot(angle2,Vstall_TO, "g", linewidth=0.8,marker = '*', label="TO")
plt.plot(angle2,Vstall_rectracted, "b", linewidth=0.8,marker = '*', label="RECTRACTED")
plt.plot(angle2,Vstall_LANDING, "r", linewidth=0.8,marker = '*', label="LANDING")
plt.xlabel('Angle in degree')
plt.ylabel('Vstall in m/s')
plt.show()




