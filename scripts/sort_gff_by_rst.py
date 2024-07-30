import os
import glob
import shutil
import datetime as dt

# get cwd and print to console
cwd = os.getcwd()
print(cwd)

sequencing_methods = ['shortread', 'longread']
rst_types = ['RST1', 'RST2', 'RST3']

# lets make sure there are directories for each sequencing method and within those, directories for each RST type
for rst_type in rst_types:
    for seq_method in sequencing_methods:
        # make the directories
        os.makedirs(f"{cwd}/{seq_method}/{rst_type}", exist_ok=True)
        print(f"{cwd}/{seq_method}/{rst_type}")

# now let's check on our annotations to make sure we've got em all! (and put em in a list for handy dandy access :)
list_of_annotations = glob.glob(pathname=f"{cwd}/*/annotation/*.gff3")

# now let's set the lists by hand to make it easier to move to remote machine :) (they're dicts but I need the keys to direct the files to the right place)
RST1 = {
    'RST1': [
          'ESI26', 'ESI26H',  'UCT109',  'UCT109H',   'UCT29',  'UCT29H',   'UCT31',
         'UCT31H',  'UCT96',  'UCT96H', 'UNY1032P',  'URI101', 'URI101H',  'URI102',
        'URI102H', 'URI107', 'URI107H',   'URI111', 'URI111H',  'URI120', 'URI120H',
          'URI39', 'URI39H',   'URI40',   'URI40H',   'URI41',  'URI41H',   'URI42',
         'URI42H',  'URI87',  'URI87H',    'URI89',  'URI89H',   'URI91',  'URI91H',
          'URI93', 'URI93H']
    }
RST2 = {
    'RST2': [
        'UCT110', 'UCT110H',  'UCT113',  'UCT113H',    'UCT30',  'UCT30H',   'UCT32',
        'UCT32H',   'UCT92',  'UCT92H', 'UNY1038P', 'UNY1090P',  'UNY149', 'UNY149P',
        'UNY208', 'UNY208P', 'UNY990P',   'URI103',  'URI103H',  'URI112', 'URI112H',
        'URI117', 'URI117H',  'URI118',  'URI118H',    'URI33',  'URI33H',   'URI44',
        'URI44H',   'URI46',  'URI46H',    'URI47',   'URI47H',   'URI48',  'URI48H',
         'URI86',  'URI86H',   'URI88',   'URI88H',   'UWI247', 'UWI247P',  'UWI248',
        'UWI248P', 'UWI263', 'UWI263P']
    }
RST3 = {
    'RST3': [
          'UCT35',  'UCT35H',  'UCT50',   'UCT50H', 'UNY1083P', 'UNY1085P', 'UNY1128P',
         'UNY169', 'UNY169P', 'UNY172',  'UNY172P',   'UNY193',  'UNY193P',   'UNY203',
        'UNY203P',   'URI34', 'URI34H',    'URI36',   'URI36H',    'URI56',   'URI56H',
         'UWI283', 'UWI283P']
    }
# pile em up
list_of_lists = [RST1, RST2, RST3]
# now let's actually move em. we'll just print the paths for now
for type in list_of_lists:
    for rst, assemblies in type.items():
        for assembly in assemblies:
            if assembly.endswith(('H', 'P')):
                # set path to longread annotation
                gff_path = f"{cwd}/longread/annotation/{assembly}.gff3"
                if os.path.isfile(gff_path):
                    print(f'{assembly}.gff3 exists. running command:')
                    print(f'os.copy(f"{cwd}/longread/annotation/{assembly}.gff3", f"{cwd}/longread/{rst}/{assembly}.gff3")')
                    # move this to the longread directory
                    shutil.copy(f"{cwd}/longread/annotation/{assembly}.gff3", f"{cwd}/longread/{rst}/{assembly}.gff3")
                else:
                    print(f'{assembly}.gff3 does not exist. skipping.')
            else:
                # set path to shortread annotation
                gff_path = f"{cwd}/shortread/annotation/{assembly}.gff3"
                if os.path.isfile(gff_path):
                    print(f'{assembly}.gff3 exists. running command:')
                    print(f'os.copy(f"{cwd}/shortread/annotation/{assembly}.gff3", f"{cwd}/shortread/{rst}/{assembly}.gff3")')
                    # move this to the shortread directory
                    shutil.copy(f"{cwd}/shortread/annotation/{assembly}.gff3", f"{cwd}/shortread/{rst}/{assembly}.gff3")
                else:
                    print(f'{assembly}.gff3 does not exist. skipping.')
print('done!')

for rst in rst_types:
    num_source = len(list_of_lists[rst_types.index(rst)][rst])
    num_destination = len(os.listdir(f"{cwd}/longread/{rst}")) + len(os.listdir(f"{cwd}/shortread/{rst}"))
    print(f'{rst} in source: {num_source}')
    print(f'{rst} in destination: {num_destination}')
    if num_source == num_destination:
        print(f'{rst} is good to go!')
    else:
        print(f'{rst}: UH OH SPAGHETTIO. something went wrong! check the output above to see what is missing.')

for rst in rst_types:
    print(f"Contents of {rst}:")
    print(os.listdir(f"{cwd}/longread/{rst}"))
    print(os.listdir(f"{cwd}/shortread/{rst}"))

print('done!')
print(dt.datetime.now())