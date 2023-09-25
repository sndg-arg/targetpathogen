import time

import pandas as pd
import requests
import sys
import json
import os
import argparse
import subprocess as sp
import tqdm
import tarfile
import shutil
import Bio.PDB as biopdb
import tempfile


class FoldSeek:
    def __init__(self):
        self.cmd_search = "curl -X POST -F q=@{filepath} -F 'mode=3diaa' -F 'database[]=pdb100' https://search.foldseek.com/api/ticket"

    def submitJob(self, filepath):
        res = json.loads(sp.check_output(self.cmd_search.format(filepath=filepath), shell=True).decode("utf-8"))
        # {"id":"PvH__CrMjPXgueMorGjjuElX2KI2IHrAcXkC2w","status":"PENDING"}
        return FoldSeekJob(res["id"])


class FoldSeekJob:
    def __init__(self, ticket, outfile=None, sleeptime=3):
        self.ticket = ticket
        self.sleeptime = sleeptime
        self.outfile = outfile
        if not self.outfile:
            self.outfile = tempfile.NamedTemporaryFile(delete=False)
            self.outfile.close()
            self.outfile = self.outfile.name
        self.df = None
        self.cmd_status = "curl -X GET https://search.foldseek.com/api/ticket/{ticket}"
        self.cmd_download = "curl -X GET https://search.foldseek.com/api/result/download/{ticket} > {outfile}"

    def wait_complete(self):
        res = json.loads(sp.check_output(self.cmd_status.format(ticket=self.ticket), shell=True).decode("utf-8"))

        # {"id":"UFMNJNBWBpGD26rdNpCfr6RorMkbzF-0GZFy1w","status":"COMPLETE"}
        if res["status"] == "ERROR":
            return Exception("error procesing FoldSeek")
        elif res["status"] == "PENDING":
            time.sleep(self.sleeptime)
            self.checkStatus()
        elif res["status"] == "COMPLETE":
            return
        else:
            return Exception("unknown FoldSeek server status:" + res["status"])

    def download_result(self):
        sp.call(self.cmd_download.format(ticket=self.ticket, outfile=self.outfile),shell=True)

    def process_result(self):
        # job.pdb	3iva_A Structure of the B12-dependent Methionine Synthase (MetH) C-teminal half with AdoHcy bound	58500	594	220	6	665	1249	1	576	1000	5,11E-59	2664	1257	576	DLAWREWPVEKRLEHALVNGITEFIEADTEEARLAAERPLHVIEGPLMAGMNVVGDLFGSGKMFLPQVVKSARVMKQAVAVLLPHMEEEKRANGGGEARESAGKILMATVKGDVHDIGKNIVGVVLACNNYEIIDLGVMVPSAKILEVAREQKVDIVGLSGLITPSLDEMAHVASELEREGFDVPLLIGGATTSRVHTAVKINPRYSLGQTVYVTDASRAVGVVSSLLSPEVRDSYKKTVRAEYLKVADAHARN---------EAEKRRLPLSQARANAFRIDWDAHQPKVPSFLGTRVFEGWDLAELARYIDWTPFFQTWELKGVFPKILDDERQGAAARQLFEDAQAMVEKIVAEAWFAPKAVIGFWPAASMGDDVRLFADEVREAELATFFTLRQQMVKRDGRPNVALADFVAPAASGKRDYVGGFVVTAGIEEVAIAERFERANDDYSSIMVKALADRFAEAFAERMHEYVRKELWGYAPDEAFTPQELIAEPYAGIRPAPGYPAQPDHTEKETLFRLLDAEAAIGVRLTESYAMWPGSSVSGLYVGHPDSYYFGVAKIERDQVEDYADRKRMSVREVERWLSPILNYVP	QAEWRSWEVNKRLEYSLVKGITEFIEQDTEEARQQATRPCEVIEGPLMDGMNVVGDLFGEGKMFLPQVVKSARVMKQAVAYLEPFIEAS------KEQCKTNGKMVIATVKGDVHDIGKNIVGVVLQCNNYEIVDLGVMVPAEKILRTAKEVNADLIGLSGLITPSLDEMVNVAKEMERQGFTIPLLIGGATTSKAHTAVKIEQNYS-GPTVYVQNASRTVGVVAALLSDTQRDDFVARTRKEYETVRIQHGRKKPRTPPVTLEA---------ARDNDFAFDWQAYTPPVAHRLGVQEVE-ASIETLRNYIDWTPFFMTWSLAGKYPRILEDEVVGVEAQRLFKDANDMLDKLSAEKTLNPRGVVGLFPANRVGDDIEIYRDETRTHVINVSHHLRQQ-TEKTGFANYCLADFVAPKLSGKADYIGAFAVTGGLEEDALADAFEAQHDDYNKIMVKALADRLAEAFAEYLHERVRKVYWGYAPNENLSNEELIRENYQGIRPAPGYPACPEHTEKATIWELLEVEKHTGMKLTESFAMWPGASVSGWYFSHPDSKYYAVAQIQRDQVEDYARRKGMSVTEVERWLAPNLGYDA	90.616,-32.713,34.188,87.804,-30.910,36.038,85.090,-33.530,36.501,82.350,-31.341,35.040,84.364,-31.163,31.818,83.822,-34.912,31.496,80.037,-34.667,31.543,77.559,-34.559,28.665,77.815,-31.352,26.610,74.302,-30.152,27.433,75.143,-31.092,31.027,78.222,-28.859,31.122,75.979,-26.057,29.869,73.514,-26.551,32.694,76.473,-26.583,35.095,78.160,-23.488,33.650,74.878,-21.577,33.640,73.966,-22.357,37.251,77.495,-21.830,38.591,77.989,-18.873,36.282,81.180,-20.189,34.673,82.068,-18.002,31.663,85.454,-19.355,30.562,84.031,-22.384,28.708,80.636,-21.153,27.487,82.081,-20.330,24.051,83.913,-23.582,23.308,81.080,-25.675,24.739,77.923,-24.367,23.096,79.435,-24.549,19.613,79.669,-28.293,20.188,76.050,-28.452,21.324,74.816,-26.685,18.195,76.368,-29.396,16.048,74.439,-32.147,17.825,71.275,-30.133,17.235,69.644,-30.049,13.806,69.150,-26.276,14.054,70.523,-23.236,15.982,67.311,-22.790,18.007,67.523,-26.216,19.632,70.485,-25.216,21.791,68.725,-22.116,23.123,65.364,-23.847,23.512,67.152,-26.958,24.754,70.050,-27.255,27.264,70.476,-23.488,27.433,67.013,-23.218,28.990,67.150,-26.476,30.958,69.842,-24.900,33.087,67.870,-21.693,33.558,64.948,-23.869,34.618,67.349,-25.299,37.197,68.320,-21.804,38.271,64.629,-21.051,38.582,63.598,-23.761,41.037,66.986,-23.652,42.728,66.231,-19.969,43.311,62.673,-20.669,44.426,63.947,-23.043,47.100,66.322,-20.562,48.697,69.547,-22.192,47.579,70.256,-19.394,45.132,69.817,-15.613,45.126,69.519,-12.665,42.741,73.231,-11.820,42.564,73.963,-15.395,41.408,70.980,-15.482,39.072,72.207,-12.497,37.053,75.548,-14.283,36.827,73.775,-17.202,35.155,72.234,-14.670,32.759,75.725,-13.488,31.870,76.529,-17.030,30.730,73.231,-17.451,28.898,73.207,-14.034,27.221,76.836,-14.679,26.256,75.952,-18.120,24.857,73.222,-16.586,22.717,75.522,-13.746,21.647,77.986,-16.228,20.156,75.113,-18.278,18.781,73.178,-15.400,17.200,75.197,-15.212,13.955,74.247,-18.847,13.311,70.551,-18.385,14.086,69.749,-15.296,12.015,72.021,-16.410,9.165,69.813,-19.511,9.174,66.471,-17.868,10.035,64.142,-18.942,12.838,63.513,-16.648,15.827,63.551,-19.487,18.389,60.559,-20.326,20.563,60.046,-17.842,23.405,57.539,-17.505,26.259,56.410,-14.070,25.150,57.676,-10.634,24.161,57.628,-7.424,26.219,58.497,-3.773,25.731,60.139,-2.260,28.801,61.056,1.389,29.322,61.555,4.118,31.883,58.829,6.738,31.566,59.484,10.040,29.743,61.526,12.507,31.782,62.646,9.905,34.309,66.188,8.549,34.271,66.701,5.869,36.938,66.588,2.322,35.549,68.326,-0.071,37.984,65.054,-1.654,39.080,63.981,-1.963,35.458,67.230,-3.653,34.417,66.561,-6.184,37.179,62.982,-6.947,36.177,64.171,-7.313,32.602,67.118,-9.649,33.207,65.019,-11.881,35.475,62.365,-12.179,32.765,64.977,-13.139,30.188,65.945,-15.880,32.638,62.460,-17.241,32.028,62.928,-17.086,28.245,60.543,-14.129,27.929,61.870,-11.709,25.300,62.523,-8.163,26.486,62.932,-5.063,24.336,64.460,-2.294,26.453,63.783,1.178,25.086,65.894,2.706,27.832,64.949,5.747,29.888,63.294,9.139,29.243,61.268,7.244,26.645,58.331,9.109,25.059,54.945,7.397,25.397,54.842,7.567,21.585,57.881,5.340,20.961,56.662,3.005,23.712,53.323,2.182,22.102,54.909,2.374,18.655,57.639,-0.139,19.515,55.157,-2.534,21.121,53.208,-2.503,17.853,56.294,-2.529,15.620,57.503,-5.813,17.125,54.256,-7.351,18.425,54.619,-7.360,22.197,52.269,-9.213,24.517,53.016,-6.850,27.373,54.182,-3.312,28.031,56.176,-2.509,31.134,56.567,0.988,32.574,59.118,2.149,35.153,59.152,5.306,37.264,61.263,6.597,40.162,59.751,10.034,40.832,56.192,11.310,41.364,56.057,13.192,38.030,56.369,9.793,36.354,52.959,8.660,37.627,51.504,11.359,35.381,53.274,10.100,32.251,52.159,6.585,33.200,48.442,7.399,33.363,49.175,9.073,30.043,50.414,5.833,28.456,47.267,4.014,29.654,44.930,6.681,28.258,46.799,6.459,24.952,46.621,2.663,24.883,42.857,3.132,25.123,42.657,5.827,22.425,44.581,3.458,20.150,42.730,0.239,20.908,45.494,-1.735,22.642,44.357,-4.986,24.246,47.582,-6.260,25.812,48.220,-6.037,29.597,49.864,-2.987,31.165,52.332,-3.611,33.986,53.332,-0.555,35.966,56.157,-0.679,38.490,58.738,1.530,40.138,59.552,3.163,43.458,56.945,5.862,42.846,54.223,3.357,41.928,51.846,1.900,44.506,48.926,-0.528,44.503,46.493,2.197,45.550,47.518,4.627,42.801,47.387,1.967,40.103,43.841,1.151,41.199,42.669,4.727,41.744,44.511,6.502,38.918,45.900,4.068,36.321,43.500,1.159,35.710,40.214,2.939,34.948,41.902,4.730,32.029,42.418,1.475,30.126,39.963,-1.203,28.967,42.785,-3.726,28.991,44.418,-5.307,32.069,46.553,-3.084,34.296,48.607,-4.769,37.000,50.925,-3.252,39.575,53.691,-5.515,40.843,55.991,-4.813,43.790,59.232,-6.433,42.623,60.946,-8.419,39.859,60.120,-11.824,41.350,56.328,-11.653,41.533,56.515,-10.137,38.056,57.483,-13.351,36.285,54.554,-14.830,38.183,52.164,-12.231,36.794,53.428,-12.339,33.219,52.944,-16.087,33.597,49.301,-16.264,34.695,48.333,-13.105,32.821,49.877,-14.524,29.640,48.443,-18.036,29.966,45.071,-19.466,28.959,43.996,-21.437,32.026,44.630,-18.521,34.385,44.644,-15.346,32.291,40.904,-14.722,32.651,40.633,-15.564,36.354,43.618,-13.359,37.243,42.261,-10.207,35.594,38.862,-10.597,37.262,40.562,-11.021,40.635,42.761,-7.994,40.005,39.705,-5.958,39.010,37.607,-6.850,42.056,40.494,-6.206,44.440,41.383,-2.970,42.648,37.797,-1.994,43.417,37.555,-2.953,47.086,40.475,-0.535,47.301,39.308,1.951,44.679,36.082,2.300,46.649,37.987,1.896,49.916,40.021,5.045,49.296,37.365,6.385,46.951,35.033,7.219,49.818,34.968,6.638,53.575,36.779,10.000,53.409,40.432,10.335,54.575,39.645,10.889,58.239,41.910,13.738,59.307,40.065,17.054,58.664,42.088,19.565,56.592,42.640,23.103,57.858,42.625,26.670,56.573,46.052,27.524,55.188,45.556,30.532,57.451,45.013,28.474,60.591,47.969,26.258,59.702,50.179,29.326,59.345,49.079,30.567,62.773,50.187,27.192,64.141,53.612,27.457,62.526,56.592,26.156,64.514,56.634,27.985,67.856,60.052,29.414,68.717,61.411,32.424,70.605,63.824,33.768,68.007,64.878,36.831,69.989,66.605,34.475,72.428,68.233,32.334,69.753,71.101,33.245,67.433,71.558,30.961,64.393,75.151,29.923,63.879,76.895,32.117,61.321,78.319,30.193,58.383,82.110,30.166,58.299,81.800,30.412,54.552,79.420,31.882,52.005,79.250,31.327,48.291,79.230,28.493,45.815,82.019,25.961,45.753,82.753,23.255,43.198,83.506,20.140,45.230,84.934,16.719,44.533,85.948,13.737,46.615,87.768,10.574,45.607,86.275,7.101,45.814,88.832,6.323,48.523,87.469,9.078,50.738,83.823,8.223,50.046,83.994,4.446,50.630,84.103,5.146,54.373,80.843,7.042,54.092,79.004,4.397,52.091,75.837,2.886,53.557,75.765,-0.667,52.169,72.319,-1.659,53.437,70.365,-0.111,50.478,72.924,-1.728,48.200,71.991,-5.200,49.519,68.288,-4.377,49.195,68.899,-3.539,45.537,70.166,-7.101,44.981,67.333,-8.764,46.909,69.515,-9.681,49.881,67.624,-8.955,53.089,69.472,-7.678,56.134,71.234,-4.657,57.575,74.988,-4.196,57.137,77.242,-5.523,58.545,75.174,-8.198,60.318,73.940,-9.402,56.909,77.437,-10.764,56.255,76.881,-13.333,59.002,73.330,-14.389,58.183,73.111,-18.175,58.033,71.155,-17.855,54.788,72.443,-14.990,52.631,75.617,-14.053,54.519,77.845,-16.035,52.176,76.293,-14.490,49.058,76.191,-11.024,50.619,79.925,-11.455,51.244,80.658,-12.374,47.642,78.466,-9.544,46.358,80.248,-7.230,48.826,83.654,-8.471,47.725,82.758,-7.910,44.072,81.487,-4.364,44.696,84.685,-3.523,46.539,86.767,-5.191,43.831,85.163,-2.846,41.309,85.629,0.181,43.544,89.322,-0.633,43.754,89.561,-0.668,39.969,87.290,2.354,39.454,89.389,4.271,41.995,92.841,3.271,40.766,91.957,4.031,37.150,89.932,7.126,38.075,87.371,5.901,35.558,84.482,7.091,37.763,84.741,10.507,39.434,82.229,11.762,42.034,81.719,15.511,41.880,79.343,18.192,43.056,78.741,21.692,44.323,77.800,23.246,47.656,76.672,26.702,48.691,75.744,28.800,51.716,73.741,32.018,51.532,72.366,34.616,53.931,68.674,33.840,54.224,65.525,34.662,56.110,61.976,33.367,56.235,59.042,35.252,54.779,55.739,33.454,55.084,56.166,29.894,53.879,59.202,30.693,51.724,62.869,31.562,52.251,64.421,34.711,50.762,68.062,34.332,49.721,70.241,37.443,49.685,73.178,38.342,47.477,75.329,39.266,50.461,75.258,40.013,54.175,73.654,43.440,53.767,70.484,41.522,52.899,69.266,44.290,50.607,69.172,42.286,47.372,67.196,39.035,47.048,68.863,36.695,44.522,66.329,33.863,44.766,63.462,32.304,46.712,63.017,28.807,48.150,59.304,28.100,48.328,58.351,25.042,50.440,54.849,23.880,51.432,53.055,23.697,54.794,51.216,20.681,56.217,47.513,21.404,55.916,46.019,19.059,58.488,44.441,19.668,61.880,46.982,19.093,64.617,46.963,19.052,68.391,50.181,19.110,70.388,52.164,19.593,67.185,52.551,22.421,64.690,51.802,22.611,60.976,55.346,22.217,59.706,56.977,24.135,56.899,60.723,23.957,56.295,60.818,27.764,56.170,59.789,27.735,59.853,63.309,26.541,60.711,64.934,29.697,59.328,64.857,32.944,61.338,61.892,35.185,60.396,62.872,38.258,58.375,60.986,40.385,60.892,63.248,39.214,63.725,66.452,40.898,62.601,68.346,37.688,63.300,70.965,36.738,60.706,70.358,33.232,59.360,71.458,31.108,56.396,70.525,28.288,54.039,72.650,25.744,52.198,72.409,23.779,48.973,74.203,20.925,47.297,74.198,18.702,44.239,76.433,15.635,43.974,76.777,12.723,41.570,78.646,9.423,41.784,80.241,10.279,38.446,79.401,11.691,35.016,80.801,8.433,33.638,78.039,6.256,35.077,76.365,5.940,31.658,79.495,5.311,29.595,80.523,2.573,32.029,77.061,1.023,32.057,77.017,0.941,28.263,80.389,-0.788,28.087,79.176,-3.678,30.249,75.885,-3.468,28.374,77.891,-4.014,25.206,79.621,-7.199,26.374,76.245,-8.480,27.626,77.443,-8.240,31.230,74.296,-7.463,33.213,76.105,-8.236,36.449,78.444,-5.244,36.621,76.186,-3.162,34.388,73.915,-3.213,37.443,76.900,-3.214,39.809,78.430,0.012,38.454,75.131,1.880,38.619,74.302,0.590,42.105,77.823,1.549,43.110,77.537,4.976,41.472,74.395,5.456,43.575,76.153,4.439,46.803,78.955,6.836,45.872,76.472,9.664,45.368,74.704,9.139,48.659,78.104,8.707,50.293,79.105,12.087,48.852,75.974,13.641,50.335,76.821,12.253,53.780,80.451,13.356,53.545,79.584,16.849,52.309,76.883,17.324,54.943,79.585,16.553,57.512,82.005,19.025,55.911,79.224,21.636,55.985,78.431,20.985,59.661,82.087,20.916,60.725,83.739,23.365,58.312,81.806,25.306,55.637,78.990,26.223,58.001,81.093,25.338,61.066,78.572,24.884,63.909,79.825,21.389,64.742,83.609,21.421,64.181,84.293,19.049,67.117,81.576,16.624,66.044,82.974,13.183,65.235,80.554,10.349,64.499,80.951,6.804,63.197,79.354,5.860,59.885,76.726,4.000,61.895,76.229,7.201,63.846,75.821,9.263,60.649,73.217,6.885,59.216,71.206,7.312,62.424,71.390,11.068,61.790,72.574,11.658,65.344,73.715,15.171,64.428,72.512,18.770,64.198,70.544,19.923,61.173,68.888,18.193,58.240,69.354,18.057,54.466,65.887,17.642,52.911,65.856,16.873,49.198,63.159,18.125,46.770,60.724,15.292,46.150,60.756,13.825,49.665,57.809,14.257,52.106,58.797,17.614,53.675,59.390,19.224,50.290,57.486,17.277,47.561,57.467,19.990,44.876,60.158,18.612,42.512,60.278,21.849,40.542,62.208,23.553,43.339,65.271,21.650,42.134,65.442,24.541,39.649,66.027,26.871,42.592,69.030,24.793,43.690,70.331,24.597,40.143,70.081,28.364,39.806,71.691,29.167,43.177,74.521,26.697,42.694,74.671,26.979,38.902,74.879,23.176,39.077,75.012,22.338,35.378,77.721,24.864,34.525,79.794,23.701,37.543,79.415,19.903,37.471,77.720,18.750,34.286,74.792,17.298,36.176,71.550,18.086,34.324,67.844,18.518,34.994,64.946,17.094,32.932,61.480,18.450,32.185,60.133,16.139,34.859,62.928,17.529,37.033,64.824,14.281,37.408,68.583,14.555,37.661,71.287,13.179,35.408,73.437,11.170,35.844,71.004,8.818,37.563,73.497,8.422,40.385,72.992,11.909,41.778,71.510,13.687,44.788,70.687,17.290,45.654,69.030,19.139,48.520,69.111,21.867,51.143,70.898,22.355,54.467,69.273,23.616,57.641,71.194,24.261,60.841,68.394,24.602,63.422,68.756,21.931,66.119,64.971,21.513,66.133,64.642,20.926,62.373,63.357,17.485,61.410,61.758,15.538,58.564,57.994,15.247,58.129,55.516,15.043,55.273,53.550,17.968,53.872,51.008,15.333,52.926,47.681,16.470,51.493,47.150,20.202,51.020,44.064,22.151,49.960,42.944,24.525,47.209,43.930,27.649,49.164,47.616,26.694,49.374,47.585,25.466,45.774,46.496,28.964,44.748,48.756,30.934,47.076,51.648,28.704,45.969,50.646,29.233,42.324,50.906,33.006,42.730,54.309,32.681,44.412,55.653,30.242,41.815,54.099,32.133,38.928,52.294,29.031,37.661,48.698,28.212,36.796,46.649,26.572,39.532,46.464,23.764,36.956,50.222,23.293,36.629,50.830,23.255,40.391,48.099,20.641,40.848,49.891,18.289,38.481,53.170,18.716,40.368,51.577,18.242,43.770,49.175,15.560,42.558,50.184,12.972,45.200,49.790,15.631,47.881,46.451,17.086,46.841,43.821,16.803,49.557,41.106,16.924,46.905,40.405,15.410,43.480,42.006,16.219,40.121	QAEWRSWEVNKRLEYSLVKGITEFIEQDTEEARQQATRPCEVIEGPLMDGMNVVGDLFGEGKMFLPQVVKSARVMKQAVAYLEPFIEASKEQCKTNGKMVIATVKGDVHDIGKNIVGVVLQCNNYEIVDLGVMVPAEKILRTAKEVNADLIGLSGLITPSLDEMVNVAKEMERQGFTIPLLIGGATTSKAHTAVKIEQNYSGPTVYVQNASRTVGVVAALLSDTQRDDFVARTRKEYETVRIQHGRKKPRTPPVTLEAARDNDFAFDWQAYTPPVAHRLGVQEVEASIETLRNYIDWTPFFMTWSLAGKYPRILEDEVVGVEAQRLFKDANDMLDKLSAEKTLNPRGVVGLFPANRVGDDIEIYRDETRTHVINVSHHLRQQTEKTGFANYCLADFVAPKLSGKADYIGAFAVTGGLEEDALADAFEAQHDDYNKIMVKALADRLAEAFAEYLHERVRKVYWGYAPNENLSNEELIRENYQGIRPAPGYPACPEHTEKATIWELLEVEKHTGMKLTESFAMWPGASVSGWYFSHPDSKYYAVAQIQRDQVEDYARRKGMSVTEVERWLAPNLGYDA	83333	Escherichia coli K-12
        cols = "x1	pdb_desc	x2	x3	x4	x5	x6	x7	x8	x9	x10	x11	x12	x13	x14	seq1    seq2    nums	seq3    tax	taxname".split()

        self.df = pd.read_csv(self.outfile, sep="\t", names=cols, index_col=False,compression='gzip')
        self.df = self.df.dropna()
        self.df["pdb"] = [x.split("_")[0] for x in self.df.pdb_desc]
        self.df["chain"] = [x.split("_")[1].split()[0] for x in self.df.pdb_desc]


    def save(self, newfile):
        shutil.copy(self.outfile, newfile)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('pdb_file', type=str)
    parser.add_argument('--wait', type=int, default=3)

    args = parser.parse_args()

    assert os.path.exists(args.pdb_file), f'"{args.pdb_file}" does not exist'

    fs = FoldSeek()
    job = fs.submitJob(args.pdb_file)
    job.wait_complete()
    job.download_result()
    job.process_result()
    sys.stderr.write(f"{job.outfile}\n")
    print(job.df)
