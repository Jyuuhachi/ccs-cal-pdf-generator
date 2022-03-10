from mbiplotter import MBIfile
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('filepath', type=str, help=r"enter a filepath ex: C:\Users\Documents\mbifile.mbi")
args = parser.parse_args()
file = args.filepath
file.encode('unicode_escape')
mbi_file = {f'{file}': {'CE07':
    {'118': {'at': 0.0, 'gauss': 0.0, 'centroid': 0.0},
    '322': {'at': 0.0, 'gauss': 0.0, 'centroid': 0.0},
    '622': {'at': 0.0, 'gauss': 0.0, 'centroid': 0.0},
    '922': {'at': 0.0, 'gauss': 0.0, 'centroid': 0.0},
    '1222': {'at': 0.0, 'gauss': 0.0, 'centroid': 0.0},
    '1522': {'at': 0.0, 'gauss': 0.0, 'centroid': 0.0},
    '1822': {'at': 0.0, 'gauss': 0.0, 'centroid': 0.0},
    '2122': {'at': 0.0, 'gauss': 0.0, 'centroid': 0.0},
    '2422': {'at': 0.0, 'gauss': 0.0, 'centroid': 0.0},
    '2722': {'at': 0.0, 'gauss': 0.0, 'centroid': 0.0},
         }
    }
}
print('fuk')
MBIfile.read_frames(MBIfile, mbi_file)
