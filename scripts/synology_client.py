#! /local/projects-t3/MSL/pipelines/packages/miniconda3/envs/interactive/bin/python3

import argparse, sys, os, glob, re, shlex
from pathlib import Path
import PyInquirer  # MUST USE PYTHON3.6 FOUND IN interactive TO USE THIS
import examples
from subprocess import run, Popen, PIPE
from enum import Enum
from datetime import datetime

class Mode(Enum):
    BROWSE = 1
    DIRECTORY_UPLOAD = 2
    DIRECTORY_DOWNLOAD = 3

SYNOLOGY_GLOBALS = {
    "START_PATH": "/volume1/NetBackup/MSL/",
    "MODE": Mode.BROWSE
}


# rsync ./directory : uploads directory and contents
# rsync ./directory/ : uploads contents
# rsync ./symlink_to_directory : uploads symlink
# rsync ./symlink_to_directory/ : uploads contents

# User sees symlink to directory: user may want to copy symlink or directory...more likely directory
# User sees symlink to file: User probably wants to copy the file...bc bash resolves the symlink before executing any commands


def main():
    user = getuser()
    while True:
        mainmenu(user) # has a choice to quit via sys.exit(0)

def getuser():
    try:
        default = os.getlogin()
    except OSError:
        default = os.environ['USER']

    questions = [
        {
            'type': 'input',
            'name': 'user_prompt',
            'message': "User: ",
            'default': default
        }
    ]
    user = PyInquirer.prompt(
        questions, style=examples.custom_style_2)['user_prompt']
    return(user)
    #lister = Popen(
    #    ["ssh", user + "@10.90.233.174", "-p", "24"],
    #    stdin=PIPE, stdout=PIPE, stderr=PIPE)

#    response = lister.stdout.readline()
 #   sys.exit(0)
  #  if re.match(".*password:", response):
   #     questions[1]['message'] = response
    #    pw = PyInquirer.prompt(
     #       questions, style=examples.custom_style_2)['pw_prompt']
      #  lister.stdin.write(pw + "\n")
       # lister.stdin.flush()

    #print(lister.stdout.readline())

   # return(lister)

def mainmenu(user):
    global SYNOLOGY_GLOBALS
    
    while True:
        questions = [
                {
                    'type': 'list',
                    'name': 'mainmenu',
                    'message': "Main menu",
                    'choices': ["DOWNLOAD", "UPLOAD", "QUIT"],
                }
            ]
        choice = inquire(questions)
        SYNOLOGY_GLOBALS["MODE"] = Mode.BROWSE # Reset state each time
        if choice == "QUIT":
            sys.exit(0)
        elif choice == "DOWNLOAD":
            download_ui(user)
        elif choice == "UPLOAD":
            choose_upload_source(user)

def choose_upload_source(user):
    wd = Path.cwd()
    global SYNOLOGY_GLOBALS
    while True:
        tree = get_local_tree(wd)
        mode_change = "DIRECTORY UPLOAD MODE" if SYNOLOGY_GLOBALS["MODE"] == Mode.BROWSE else "BROWSE MODE"
        tree = ["CANCEL", mode_change] + tree[1:len(tree)]
        questions = [
            {
                'type': 'list',
                'name': 'chooser',
                'message': "Please choose upload source. MODE: " + SYNOLOGY_GLOBALS["MODE"].name,
                'choices': tree,
            }
        ]
        choice = PyInquirer.prompt(questions, style=examples.custom_style_2)['chooser']
        
        def process_directory(path):
            nonlocal wd
            if SYNOLOGY_GLOBALS["MODE"] == Mode.BROWSE:
                wd = path

            if SYNOLOGY_GLOBALS["MODE"] == Mode.DIRECTORY_UPLOAD:
                try:
                    choose_upload_dest(user, path, "directory")
                except ConnectionError as e:
                    print(e)
                    return
                except Exception as e:
                    print(e + "\n")
        # FILE OR LINK
        if choice[0] == "-":
            basename = choice.split()[8]
            fullname = wd / Path(basename)
            choose_upload_dest(user, fullname, "file")

        if choice[0] == "l":
            target = Path(choice.split()[8]).resolve()
            if target.is_dir():
                process_directory(target)
            else:
                choose_upload_dest(user, target, "file")
                    
        # GO BACK TO MAIN MENU
        elif choice == "CANCEL":
            break

        # DIRECTORY -> navigate to that directory, with . and .. resolution
        elif choice[0] == 'd':
            basename = choice.split()[8]
            fullname = (wd / Path(basename)).resolve()
            process_directory(fullname)
                    
        elif choice == "DIRECTORY UPLOAD MODE":
            SYNOLOGY_GLOBALS["MODE"] = Mode.DIRECTORY_UPLOAD
        elif choice == "BROWSE MODE":
            SYNOLOGY_GLOBALS["MODE"] = Mode.BROWSE
            
    return

def inquire(questions):
    answers = PyInquirer.prompt(
            questions,
            style=examples.custom_style_2)
    for answer in answers.values():
        return answer # only one answer expected

# may throw ConnectionError
def choose_upload_dest(user, source, ul_type):
    global SYNOLOGY_GLOBALS
    wd = Path(SYNOLOGY_GLOBALS["START_PATH"])
    success = True 

    while True:
        tree = get_remote_tree(user, wd) # may throw ConnectionRefusedError
        # choices.sort()
        tree = ["CANCEL", "UPLOAD HERE"] + tree[1:len(tree)]
        questions = [
            {
                'type': 'list',
                'name': 'chooser',
                'message': "Please choose upload destination:",
                'choices': tree,
            }
        ]
        choice = PyInquirer.prompt(questions, style=examples.custom_style_2)['chooser']

        # FILE TO FILE: Upload and overwrite
        if choice[0] == "-" and ul_type == "file":
            basename = choice.split()[8]
            fullname = wd / Path(basename)
            if(type(source) == str):
                if ask_upload(user, source, fullname, "file to file"):
                    break
        
        # Can't upload directory to a file
        elif choice[0] == "-" and ul_type == "directory":
            raise Exception("Can't upload a directory to a file\n")

        elif choice[0] == "l":
            target = Path(choice.split()[8])
            if target.is_dir():
                wd = (wd / Path(target)).resolve()
            else:
                raise Exception("Can't upload to a symlink\n")
        
        # GO BACK TO MAIN MENU
        elif choice == "CANCEL":
            success = False
            break

        # DIRECTORY: resolve absolute path, and figure out whether to enter or upload
        elif choice[0] == 'd':
            basename = choice.split()[8]
            fullpath = (wd / Path(basename)).resolve()
            wd = fullpath

        elif choice == "UPLOAD HERE":
            if ask_upload(user, list(source), wd, "%s to directory" % ul_type):
                break

        elif choice == "BROWSE MODE":
            SYNOLOGY_GLOBALS["MODE"] = Mode.BROWSE
    return success

def archive_project(user, name, subdirs):
    # make target directory
    dest = Path("/") / Path("volume1") / Path("NetBackup") / Path("MSL") / Path("projects")
    year = str(datetime.today().year)
    dest = dest / Path(f"{year}_{name.upper()}")
    cmd = ["ssh", user + "@10.90.233.174", "-o", "ConnectTimeout=10", "-p", "24", "mkdir", "-p", str(dest)]
    print("Creating project directory on Rosalind...")
    print(" ".join(cmd))
    # process = Popen(
    #         cmd,
    #         stdin=PIPE, stdout=PIPE, stderr=PIPE)
    # stderr = process.stderr.read().decode(encoding='utf8')
    # print(stderr)
    # if "Connection timed out" in stderr:
    #     raise ConnectionError(stderr)
    # elif len(stderr) > 0:
    #     raise Exception(stderr)

    # upload all
    full_paths = [str(Path(subdir).resolve()) for subdir in subdirs]
    for dir in full_paths:
        if not Path(dir).exists():
            raise FileNotFoundError(dir + " does not exist")
        if not Path(dir).is_dir():
            raise Exception(dir + " is not a directory")
    ask_upload(user, full_paths, dest, "directory to directory")

def ask_upload(user, source, dest, ul_type):
    if ul_type == "directory to directory":
        msg = "Upload the following directory/ies\n\n%s\n\nto %s now?"
    elif ul_type == "file to directory":
        msg = "Upload the following file(s)\n\n%s\n\nto %s now?"
    elif ul_type == "file to file":
        msg = "Overwrite file %s to %s now?"

    msg = msg % ("\n".join(source), dest)
    questions = [
        {
            'type': 'confirm',
            'name': 'confirm_file',
            'message': msg,
            'default': False
        }
    ]
    ok = inquire(questions)

    if ok:
        cmd = "rsync -azh --append-verify --progress %s %s@10.90.233.174:%s" % (" ".join(source), user, dest)
        print(cmd)
        run(shlex.split(cmd))
    return ok
    
# path is automatically cast to str
# may throw ConnectionError
def get_remote_tree(user, path=SYNOLOGY_GLOBALS["START_PATH"]):
    cmd = ["ssh", user + "@10.90.233.174", "-o", "ConnectTimeout=10", "-p", "24", "ls", "-al", "--group-directories-first", str(path)]
    print(" ".join(cmd))
    process = Popen(
            cmd,
            stdin=PIPE, stdout=PIPE, stderr=PIPE)
    stdout = process.stdout.read().decode(encoding='utf8').split("\n")
    stderr = process.stderr.read().decode(encoding='utf8')
    print(stderr)
    if "Connection timed out" in stderr:
        raise ConnectionError(stderr)
    elif len(stderr) > 0:
        raise Exception(stderr)
    else:
        #return 
        return stdout[0:len(stdout) - 1] # last line of ls stdout is blank

def get_local_tree(path):
    cmd = ["ls", "-al", "--group-directories-first", str(path)]
    print(" ".join(cmd))
    process = Popen(cmd, stdin=PIPE, stdout=PIPE, stderr=PIPE)
    lines = process.stdout.read().decode(encoding='utf8').split("\n")

    return lines[0:len(lines) - 1]# last line of ls stdout is blank

def download_ui(user):
    global SYNOLOGY_GLOBALS
    wd = Path(SYNOLOGY_GLOBALS["START_PATH"])
    while True:
        try:
            tree = get_remote_tree(user, wd)
        except ConnectionError as e:
            print(e)
            return
        # choices.sort()
        mode_change = "DIRECTORY DOWNLOAD MODE" if SYNOLOGY_GLOBALS["MODE"] == Mode.BROWSE else "BROWSE MODE"
        tree = ["CANCEL", mode_change] + tree[1:len(tree)]
        questions = [
            {
                'type': 'list',
                'name': 'chooser',
                'message': "MODE: " + SYNOLOGY_GLOBALS["MODE"].name,
                'choices': tree,
            }
        ]
        choice = PyInquirer.prompt(questions, style=examples.custom_style_2)['chooser']

        def process_directory(fullpath):
            nonlocal wd
            if SYNOLOGY_GLOBALS["MODE"] == Mode.BROWSE:
                wd = fullpath

            if SYNOLOGY_GLOBALS["MODE"] == Mode.DIRECTORY_DOWNLOAD:
                ask_download(user, fullpath, "directory")

        # FILE
        if choice[0] == "-":
            basename = choice.split()[8]
            fullname = wd / Path(basename)
            ask_download(user, fullname, "file")
            
        if choice[0] == "l":
            # Determine if symlink is to directory or file
            target = Path(choice.split()[8]).resolve()
            if target.is_dir():
                process_directory(target)
            else:
                ask_download(user, target, "file")

        # GO BACK TO MAIN MENU
        elif choice == "CANCEL":
            break

        # DIRECTORY -> navigate to that directory, with . and .. resolution
        elif choice[0] == 'd':
            basename = choice.split()[8]
            fullpath = (wd / Path(basename)).resolve()
            process_directory(fullpath)

        elif choice == "DIRECTORY DOWNLOAD MODE":
            SYNOLOGY_GLOBALS["MODE"] = Mode.DIRECTORY_DOWNLOAD
        elif choice == "BROWSE MODE":
            SYNOLOGY_GLOBALS["MODE"] = Mode.BROWSE
                    
    return

def ask_download(user, filepath, dl_type):
    
    questions = [
                {
                    'type': 'confirm',
                    'name': 'confirm_file',
                    'message': 'Download ' + dl_type + " " + str(filepath.name) + ' now?',
                    'default': False
                }
            ]
    ok = PyInquirer.prompt(
        questions,
        style=examples.custom_style_2)['confirm_file']
    if ok:
        cmd = "rsync -azh --append-verify --progress " + user + "@10.90.233.174:" + str(filepath) + " ./"
        print(cmd)
        run(shlex.split(cmd))
    return

########################################################################################
# ARGUMENTS
########################################################################################

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""Interactive console for transferring runs to and from Synology.

        To make it work without a password, please follow these instructions:
            
    Do this on thanos or magog: https://help.github.com/articles/generating-a-new-ssh-key-and-adding-it-to-the-ssh-agent/#generating-a-new-ssh-key
    Then:
    rsync -a ~/.ssh/id_rsa.pub SYNOLOGY_USER@10.90.233.174:/var/services/homes/SYNOLOGY_USER

    ssh SYNOLOGY_USER@10.90.233.174 -p 24
    cd /var/services/homes/
    chmod 755 SYNOLOGY_USER
    cd SYNOLOGY_USER
    mkdir .ssh
    chmod 700 .ssh
    touch ~/.ssh/authorized_keys
    cat ~/id_rsa.pub >> ~/.ssh/authorized_keys
    rm ~/id_rsa.pub
    chmod 600 authorized_keys

    Ensure this command works:
    rsync --list-only SYNOLOGY_USER@10.90.233.174:/volume1/NetBackup/MSL


            """
    )


    # parser.add_argument('-p', metavar="project name", help="Project name (REQUIRED in coassembly and merged modes). User will be prompted to delete any existing directory with this name.")
    # parser.add_argument('-s', '--samples', help="Name of samples file to create", default = "SQM_manifest.txt")
    # parser.add_argument('--exclude', nargs="*", help="Directories that do not contain sample raw reads.", default=[])
    # parser.add_argument('--read-dir', metavar="DIR", help="Change mode to find all read files in <DIR>, named <sampleA>_1.{fastq,fasta}[.gz], <sampleA>_2.{fastq,fasta}[.gz], etc. ALSO overrides any present -f|-seq for convenience.")
    # # parser.add_argument('--no-review', action='store_true', help="Execute SqueezeMeta without manually reviewing sample file.")
    # parser.add_argument('--reads-prefix', default='', help="Only find raw files with this prefix")

    args = parser.parse_args()
    
    main() 
