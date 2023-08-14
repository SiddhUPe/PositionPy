"""All formulas from https://aa.quae.nl/en/reken/hemelpositie.html and https://aa.quae.nl/en/reken/zonpositie.html"""
from PositionPy.julian import *
from astroquery.simbad import Simbad
from astroquery.exoplanet_orbit_database import ExoplanetOrbitDatabase
from Solarflare.solarflare import *
from Solarflare.nightcalc import *


def sinD(deg): return Decimal(math.sin(math.radians(deg)))
def tanD(deg): return Decimal(math.tan(math.radians(deg)))
def cosD(deg): return Decimal(math.cos(math.radians(deg)))
def asinD(deg): return Decimal(math.asin(math.radians(deg)))
def atanD(deg): return Decimal(math.atan(math.radians(deg)))
def acosD(deg): return Decimal(math.acos(math.radians(deg)))
def atanD2(deg, deg2): return Decimal(math.atan2(math.radians(deg), math.radians(deg2)))


class invalidPlanetError(Exception):
    def __init__(self, message="Planet is not recognized. Check if it is one of 9 ones in our Solar System"):
        self.message = message
        super().__init__(self.message)

class notFoundInDatabaseError(Exception):
    def __init__(self, message="Not found in the SIMBAD database - please try something else."):
        self.message = message
        super().__init__(self.message)

def mean_anomaly(planet="earth", date=datetime.datetime.now()):
    planetMCONST = {
        "mercury": (174.7948, 4.09233445),
        "venus": (50.4161, 1.60213034),
        "earth": (357.5291, 0.98560028),
        "mars": (19.3730, 0.52402068),
        "jupiter": (20.0202, 0.08308529),
        "saturn": (317.0207, 0.03344414),
        "uranus": (141.0498, 0.01172834),
        "neptune": (256.2250, 0.00598103),
        "pluto": (14.882, 0.00396),
    }
    if not planetMCONST.get(planet, None):
        # incorrect planet entered
        raise invalidPlanetError

    else:
        M0, M1 = planetMCONST[planet]
        decM0 = Decimal(M0)
        decM1 = Decimal(M1)
        J = julian_date(date)
        J2000 = 2451545
        M = (decM0 + decM1 * (J - Decimal(J2000))) % Decimal(360)
        return M

def ecliptical_longitude(planet="earth", date=datetime.datetime.now()):
    planetPERIHELION = {
        "mercury": 230.3265,
        "venus": 73.7576,
        "earth": 102.9373,
        "mars": 71.0041,
        "jupiter": 237.1015,
        "saturn": 99.4587,
        "uranus": 5.4634,
        "neptune": 182.2100,
        "pluto": 184.5484,
    }

    if not planetPERIHELION.get(planet, None):
        raise invalidPlanetError

    else:
        M = mean_anomaly(planet, date)
        II = planetPERIHELION[planet]
        C = equation_of_center(planet, date)
        long = Decimal(M) + Decimal(II) + C + Decimal(180)
        return long % Decimal(360)


def equation_of_center(planet="earth", date=datetime.datetime.now()):
    planetCCONST = {
        "mercury": [23.4400, 2.9818, 0.5255, 0.1058, 0.0241, 0.0055],  # 0.0026 is the maximum error
        "venus": [0.7758, 0.0033, 0, 0, 0, 0],  # 0.0000 is the maximum error
        "earth": [1.9148, 0.0200, 0.0003, 0, 0, 0],  # 0.0000 is the maximum error
        "mars": [10.6912, 0.6228, 0.0503, 0.0046, 0.0005, 0],  # 0.0001 is the maximum error
        "jupiter": [5.5549, 0.1683, 0.0071, 0.0003, 0, 0],  # 0.0001 is the maximum error
        "saturn": [6.3585, 0.2204, 0.0106, 0.0006, 0, 0],  # 0.0001 is the maximum error
        "uranus": [5.3042, 0.1534, 0.0062, 0.0003, 0, 0],  # 0.0001 is the maximum error
        "neptune": [1.0302, 0.0058, 0, 0, 0, 0],  # 0.0001 is the maximum error
        "pluto": [28.3150, 4.3408, 0.9214, 0.2235, 0.0627, 0.0174]  # 0.0096 is the maximum error
    }

    if not planetCCONST.get(planet, None):
        # incorrect planet entered
        raise invalidPlanetError

    else:
        c_row = planetCCONST[planet]
        c1 = Decimal(c_row[0])
        c2 = Decimal(c_row[1])
        c3 = Decimal(c_row[2])
        c4 = Decimal(c_row[3])
        c5 = Decimal(c_row[4])
        c6 = Decimal(c_row[5])
        M = mean_anomaly(planet, date)
        C = c1 * Decimal(sinD(M)) + c2 * Decimal(sinD(2 * M)) + \
            c3 * Decimal(sinD(3 * M)) + c4 * Decimal(sinD(4 * M)) + \
            c5 * Decimal(sinD(5 * M)) + c6 * Decimal(sinD(6 * M))
        return C

def true_anomaly(planet="earth", date=datetime.datetime.now()):
    C = equation_of_center(planet, date)
    M = mean_anomaly(planet, date)
    return C + M

def sun_distance(planet="earth", date=datetime.datetime.now()):
    planetAUCONST = {
        "mercury": 0.37073,
        "venus": 0.72330,
        "earth": 0.99972,
        "mars": 1.51039,
        "jupiter": 5.19037,
        "saturn": 9.52547,
        "uranus": 19.17725,
        "neptune": 30.10796,
        "pluto": 37.09129,
    }

    if not planetAUCONST.get(planet, None):
        # incorrect planet entered
        raise invalidPlanetError

    else:
        planetECONST = {
            "mercury": 0.20563,
            "venus": 0.00677,
            "earth": 0.01671,
            "mars": 0.09340,
            "jupiter": 0.04849,
            "saturn": 0.05551,
            "uranus": 0.04630,
            "neptune": 0.00899,
            "pluto": 0.2490,
        }
        e = Decimal(planetECONST.get(planet))
        a = Decimal(planetAUCONST.get(planet))
        v = Decimal(true_anomaly(planet, date))
        r = a/(Decimal(1) + e * Decimal(cosD(v)))
        return r


def helio_ecliptical_coordinates(planet="earth", date=datetime.datetime.now()):
    planetOMEGACONST = {
        "mercury": 48.331,
        "venus": 76.680,
        "earth": 174.873,
        "mars": 49.558,
        "jupiter": 100.464,
        "saturn": 113.666,
        "uranus": 74.006,
        "neptune": 131.784,
        "pluto": 110.307,
    }
    planetSOMEGACONST = {
        "mercury": 29.125,
        "venus": 54.884,
        "earth": 288.064,
        "mars": 286.502,
        "jupiter": 273.867,
        "saturn": 339.391,
        "uranus": 98.999,
        "neptune": 276.340,
        "pluto": 113.768,
    }
    planetICONST = {
        "mercury": 7.005,
        "venus": 3.395,
        "earth": 0.000,
        "mars": 1.850,
        "jupiter": 1.303,
        "saturn": 2.489,
        "uranus": 0.773,
        "neptune": 1.770,
        "pluto": 17.140,
    }

    w = Decimal(planetSOMEGACONST.get(planet, Decimal(288.064)))
    W = Decimal(planetOMEGACONST.get(planet, Decimal(174.873)))
    i = planetICONST.get(planet, Decimal(0))
    r = Decimal(sun_distance(planet, date))
    v = Decimal(true_anomaly(planet, date))

    x = r * ((cosD(W) * cosD(w + v)) - sinD(W) * cosD(i) * sinD(w + v))
    y = r * ((sinD(W) * cosD(w + v)) + cosD(W) * cosD(i) * sinD(w + v))
    z = r * sinD(i) * sinD(w + v)

    return Decimal(x), Decimal(y), Decimal(z)

def geo_ecliptical_coordinates(planet="jupiter", date=datetime.datetime.now()):
    xp, yp, zp = helio_ecliptical_coordinates(planet, date)
    xe, ye, ze = helio_ecliptical_coordinates("earth", date)
    x = xp - xe
    y = yp - ye
    z = zp - ze
    return x, y, z

def query_simbad(obj_name):
    """
    :param obj_name: name of the SIMBAD object
    :return: the RA and declination of the object
    """
    try:
        # Query the SIMBAD database for the object
        result_table = Simbad.query_object(obj_name)
        # Extract the RA and Dec coordinates from the result table
        ra = result_table['RA'][0]
        dec = result_table['DEC'][0]
        # Return the RA and Dec coordinates as floats
        return ra, dec
    except:
        raise notFoundInDatabaseError

# 103 common stars
def polaris(): return query_simbad("polaris")
def sirius(): return query_simbad("sirius")
def betelgeuse(): return query_simbad("betelgeuse")
def vega(): return query_simbad("vega")
def arcturus(): return query_simbad("arcturus")
def capella(): return query_simbad("capella")
def alphard(): return query_simbad("alphard")
def eta_leonis(): return query_simbad("eta leonis")
def regulus(): return query_simbad("regulus")
def eta_draconis(): return query_simbad("eta draconis")
def alderamin(): return query_simbad("alderamin")
def scheat(): return query_simbad("scheat")
def alpheratz(): return query_simbad("alpheratz")
def almach(): return query_simbad("almach")
def markab(): return query_simbad("markab")
def fomalhaut(): return query_simbad("fomalhaut")
def alnair(): return query_simbad("alnair")
def gamma_gruis(): return query_simbad("gamma gruis")
def peacock(): return query_simbad("peacock")
def shaula(): return query_simbad("shaula")
def nunki(): return query_simbad("nunki")
def altair(): return query_simbad("altair")
def rasalhague(): return query_simbad("rasalhague")
def sabik(): return query_simbad("sabik")
def antares(): return query_simbad("antares")
def kornephoros(): return query_simbad("kornephoros")
def alkaid(): return query_simbad("alkaid")
def alioth(): return query_simbad("alioth")
def dubhe(): return query_simbad("dubhe")
def kochab(): return query_simbad("kochab")
def gamma_cassiopeiae(): return query_simbad("gamma cassiopeiae")
def mirach(): return query_simbad("mirach")
def algenib(): return query_simbad("algenib")
def canopus(): return query_simbad("canopus")
def alpha_centauri(): return query_simbad("alpha centauri")
def pi_capricorni(): return query_simbad("pi capricorni")
def ascella(): return query_simbad("ascella")
def eta_carinae(): return query_simbad("eta carinae")
def eltanin(): return query_simbad("eltanin")
def cebalrai(): return query_simbad("cebalrai")
def zeta_ophiuchi(): return query_simbad("zeta ophiuchi")
def zeta_centauri(): return query_simbad("zeta centauri")
def eta_centauri(): return query_simbad("eta centauri")
def vindemiatrix(): return query_simbad("vindemiatrix")
def delta_cygni(): return query_simbad("delta cygni")
def porrima(): return query_simbad("porrima")
def denebola(): return query_simbad("denebola")
def mizar(): return query_simbad("mizar")
def megrez(): return query_simbad("megrez")
def phecda(): return query_simbad("phecda")
def chi_ursae_majoris(): return query_simbad("chi ursae majoris")
def psi_ursae_majoris(): return query_simbad("psi ursae majoris")
def tania_australis(): return query_simbad("tania australis")
def tania_borealis(): return query_simbad("tania borealis")
def theta_ursae_majoris(): return query_simbad("theta ursae majoris")
def kappa_ursae_majoris(): return query_simbad("kappa ursae majoris")
def iota_ursae_majoris_a(): return query_simbad("iota ursae majoris a")
def muscida(): return query_simbad("muscida")
def upsilon_ursae_majoris(): return query_simbad("upsilon ursae majoris")
def beta_leonis_minoris(): return query_simbad("beta leonis minoris")
def kappa_scorpii(): return query_simbad("kappa_scorpii")
def sargas(): return query_simbad("sargas")
def epsilon_scorpii(): return query_simbad("epsilon scorpii")
def n_scorpii(): return query_simbad("n scorpii")
def h_scorpii(): return query_simbad("h scorpii")
def k_scorpii(): return query_simbad("k scorpii")
def zeta_scorpii(): return query_simbad("zeta scorpii")
def eta_scorpii(): return query_simbad("eta scorpii")
def polis(): return query_simbad("polis")
def phi_sagittarii(): return query_simbad("phi sagittarii")
def tau_sagittarii(): return query_simbad("tau sagittarii")
def iota_sagittarii(): return query_simbad("iota sagittarii")
def rukbat(): return query_simbad("rukbat")
def mu_coronae_australis(): return query_simbad("mu coronae australis")
def alpha_telescopii(): return query_simbad("alpha telescopii")
def zeta_telescopii(): return query_simbad("zeta telescopii")
def delta_telescopii(): return query_simbad("delta telescopii")
def theta_arae(): return query_simbad("theta arae")
def alpha_arae(): return query_simbad("alpha arae")
def beta_arae(): return query_simbad("beta arae")
def eta_arae(): return query_simbad("eta arae")
def hamal(): return query_simbad("hamal")
def sheratan(): return query_simbad("sheratan")
def alpha_trianguli(): return query_simbad("alpha trianguli")
def beta_trianguli(): return query_simbad("beta trianguli")
def gamma_trianguli(): return query_simbad("gamma trianguli")
def gamma_arietes(): return query_simbad("gamma arietes")
def menkar(): return query_simbad("menkar")
def hercules_109(): return query_simbad("109 hercules")
def HD169191(): return query_simbad("HD 169191")
def HD168270(): return query_simbad("HD 168270")
def HD169222(): return query_simbad("HD 169222")
def HD167008(): return query_simbad("HD 167008")
def iq_hercules(): return query_simbad("iq hercules")
def vulpeculae1(): return query_simbad("1 Vulpeculae")
def herculis93(): return query_simbad("93 Herculis")
def HD166095(): return query_simbad("HD 166095")
def ophiuchi72(): return query_simbad("72 Ophiuchi")
def alpha_herculis(): return query_simbad("alpha herculis")
def zeta_herculis(): return query_simbad("zeta herculis")
def eta_herculis(): return query_simbad("eta herculis")
def pi_herculis(): return query_simbad("pi herculis")
def rho_herculis(): return query_simbad("rho herculis")
def c_herculis(): return query_simbad("c herculis")


# 10 exoplanets that are very similar to earth
def gliese667Cc(): return query_simbad("Gliese 667Cc")
def kepler22b(): return query_simbad("Kepler-22b")
def kepler69c(): return query_simbad("Kepler-69c")
def kepler62f(): return query_simbad("Kepler-62f")
def kepler186f(): return query_simbad("Kepler-186f")
def kepler442b(): return query_simbad("Kepler-442b")
def kepler452b(): return query_simbad("Kepler-452b")
def kepler1649c(): return query_simbad("Kepler-1649c")
def proxima_centauri_b(): return query_simbad("Proxima Centauri b")
def TRAPPIST1e(): return query_simbad("TRAPPIST-1e")

# 5 Nebulas that are commonly known
def helix_nebula(): return query_simbad("Helix Nebula")
def orion_nebula(): return query_simbad("Orion Nebula")
def crab_nebula(): return query_simbad("Crab Nebula")
def horsehead_nebula(): return query_simbad("Horsehead Nebula")
def carina_nebula(): return query_simbad("Carina Nebula")

# alternate method of finding the declination and right ascension of some exoplanets
def exoplanet_coordinates(exoplanet):
    planetTable = ExoplanetOrbitDatabase.query_planet(exoplanet)
    return planetTable['sky_coord']

class Astronomer:
    # I am using Solarflare, another python library I have programmed

    def __init__(self, latitude, longitude, height=0):
        self.lat = latitude
        self.long = longitude
        self.altitude = height

    @staticmethod
    def solarDeclination(date=datetime.datetime.now()):
        return declination(date)

    @staticmethod
    def solarRightAscension(date=datetime.datetime.now()):
        return right_ascension(date)

    def sunriseTime(self, date=datetime.datetime.now()):
        return sunrise(self.lat, self.long, date)

    def sunsetTime(self, date=datetime.datetime.now()):
        return sunset(self.lat, self.long, date)

    def hourAngle(self, date=datetime.datetime.now()):
        return hour_angle(self.long, date)

    def solarAltitude(self, date=datetime.datetime.now()):
        return altitude(self.lat, self.long, date)

    def solarAzimuth(self, date=datetime.datetime.now()):
        return azimuth(self.lat, self.long, date)

    @staticmethod
    def lunarDeclination(date=datetime.datetime.now()):
        return lunar_dec(date)

    @ staticmethod
    def lunarRightAscension(date=datetime.datetime.now()):
        return lunar_ra(date)

    @staticmethod
    def lunarPhaseInformation(date=datetime.datetime.now()):
        information = {
            "phase": get_lunar_phase(date),
            "age": get_lunar_age(date),
            "age_percent": get_lunar_age_percent(date),
            "latitude": lunar_coordinates()[0],
            "longitude": lunar_coordinates()[1],
            "distance": lunar_coordinates()[2],
        }
        return information
