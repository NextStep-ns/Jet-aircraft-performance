
Max_weight = 3920
Wing_area = 19.5
Max_Lift_Coeff_Flap_Retracted = 1.52 # for alpha equal to 20 degrees
Max_Lift_Coeff_Flap_Partially_Extended_TO = 1.93
Max_Lift_Coeff_Flap_Partially_Extended_LANDING = 2.2
Maximum_Thrust_Of_Each_Turbojet = 6490
Fuel_Capacity = 1110
Specifi_Consumption = 0.024

def Feet_Into_Meter(length):
    return length*0.3048   #length in feet initially

def Meter_Into_Feet(length):
    return length/0.3048

def Lift(alpha):
    return 0.02+0.075*alpha

def Polar(lift):
    return 0.025+(0.05*(lift**2))

def Pressure(altitude):
    return 101325*(1-(22.557*(10**-6)*altitude))**5.225

def Temperature(altitude):
    return 288-(6.5*altitude*(10**-3))

################## CLIMBING ###############
r = 0.1 # coeff of friction
V_lof = 46 #flight speed in m/s
Take_Off_Height = 1900 # Expressed in feet
g = 9.81

############## CLIMBING-Q1 ###################
v_climbing = 88
vario = Feet_Into_Meter(3000)/60

############## CRUISING FLIGHT -Q1 ##############
epsilon = 0.05
Cd0 = 0.025
gamma_lifting_ceiling = 1.4
mach_number = 0.34

############## CRUISING FLIGHT -Q4 ##############
cs = (0.093/3600)