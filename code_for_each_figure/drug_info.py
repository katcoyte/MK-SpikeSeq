# Author: Katharine Z. Coyte
#
# License: Apache License 2.0
#
# Sort any drugs taken into antibacterials / antifungals / vaccines


import pandas as pd
import numpy as np

def load_drug_types():
      antibacterials = ['Amoxicillin',
                        'Ampicillin',
                        'Ampicillin/Sulbactam (Unasyn)',
                        'Bacitracin',
                        'Cefazolin',
                        'Cefotaxime',
                        'Ceftazidime',
                        'Ceftriaxone',
                        'Cephalexin',
                        'Clindamycin',
                        'Erythromycin',
                        'Erythromycin ophthalmic',
                        'Gentamicin',
                        'Gentamicin ophthalmic',
                        'Meropenem',
                        'Mupirocin (Bactroban)',
                        'Oxacillin',
                        'Penicillin G',
                        'Rifampin',
                        'Tobramycin/Dexamethasone',
                        'Vancomycin',
                        'Zosyn (Piperacillin/tazobactam)']

      antifungals = ['Amphotericin B',
               'Amphotericin B liposomal',
               'Fluconazole',
               'Miconazole Powder',
               'Nystatin',
               'Nystatin Ointment',]

      vaccines = ['DTaP Vaccine (120 day)',
            'DTaP Vaccine (60 day)',
            'Haemophilus b (Hib) vaccine (60 day)',
            'Hepatitis B vaccine',
            'Hepatitis B vaccine (120 day)',
            'Hepatitis B vaccine (60 day)',
            'IP vaccine (Inactivated polio) (120 day)',
            'IP vaccine (Inactivated polio) (60 day)',
            'Pentacel (DTaP/Hib/IPV combination vaccine) (120 day)',
            'Pentacel (DTaP/Hib/IPV combination vaccine) (60 day)',
            'Pneumococcal 13 vaccine (Prevnar) (120 day)',
            'Pneumococcal 13 vaccine (Prevnar) (60 day)']

      return(antibacterials, antifungals, vaccines)
