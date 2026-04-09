import shutil
from astropy.time import Time

# path from NRAO archive with messy pipeline numbers etc.
mspath = "/Users/jimmylynch/Desktop/radio/observations/25A-060/pipeline.61139.014699073974/25A-060.sb48980598.eb49031908.60878.251126539355.ms"

# This function turns your NRAO formatted mspath (eg, "program/pipeline.61139.014699073974/25A-060.sb48980598.eb49031908.60878.251126539355.ms") to
# "program/target/program.target.date/program.target.date.ms" which is cleaner and used by the script. If your .ms is already in that format, it will
# throw an error to not overwrite it, so only run this once!
def rename_ms(mspath, field=2):

    # old directories
    pipeline_dir = os.path.dirname(mspath)
    program_dir = os.path.dirname(pipeline_dir)

    # listobs
    lobs = listobs(mspath)    

    # scrape and format start time
    obs_date = lobs["BeginTime"]
    t = Time(obs_date, format='mjd')
    obs_date_formatted = t.strftime('%Y-%m-%d')

    # get program name from pipeline path
    program = os.path.basename(program_dir)

    # scrape target: assuming field = 2 but can change in function inputs
    target = lobs[f"field_{field}"]["name"]  

    # new directories: going for /program/target/program.target.date/program.target.date.ms
    source_dir = f"{program_dir}/{target}"
    os.makedirs(source_dir, exist_ok=True)
    obs_dir = f"{source_dir}/{program}.{target}.{obs_date_formatted}"
    os.makedirs(obs_dir, exist_ok=False)
    
    new_mspath = f"{obs_dir}/{program}.{target}.{obs_date_formatted}.ms"
    shutil.move(mspath, new_mspath)
    print(f"Moved ms to {new_mspath}")

    return target, program, obs_date_formatted

rename_ms(mspath)