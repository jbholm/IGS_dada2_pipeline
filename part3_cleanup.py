#! /local/projects-t3/MSL/pipelines/packages/miniconda3/envs/interactive/bin/python3

import argparse, sys, os, glob, re, shlex, json, shutil, gzip, getpass, stat, tempfile, time
from pathlib import Path # Use of Paths as os.pathlike objects requires Python3.6
import PyInquirer  # MUST USE PYTHON3.6 FOUND IN interactive TO USE THIS
import examples
from subprocess import run, Popen, PIPE
from enum import Enum
import scripts.synology_client as synology  # MUST BE IN SAME DIRECTORY AS THIS SCRIPT
from jira import JIRA
from jira.exceptions import JIRAError

pipeline_dir = os.path.dirname(os.path.realpath(__file__))

with (Path(pipeline_dir) / Path("config.json")).open("r") as fh:
    CONFIG = json.load(fh)

pathfinder = {
    "intermediate_taxonomies": Path("TAXONOMY_INTERMEDIATES"),
    "final_taxonomies": Path("TAXONOMY_FINAL"),
    "counts-by-taxon": Path("COUNTS_BY_TAXON"),
    "counts-by-ASV": Path("COUNTS_BY_ASV")
}


def main(args):
    args.project = choose_project(default=args.project)

    while True:
        questions = [
            {
                "type": "list",
                "name": "chooser",
                "message": "Choose operation:",
                "choices": [
                    {
                        "name": "PREPARE FOR DELIVERY - Reorganize files and delete pipeline artifacts",
                        "value": "D",
                    },
                    {
                        "name": "UPLOAD TO JIRA - Upload select files to MSL JIRA ticket",
                        "value": "J",
                    },
                    # {
                    #     'name': "ARCHIVE - copy to Rosalind",
                    #     'value': 'A'
                    # },
                    {"name": "BACK", "value": "B"},
                    {"name": "QUIT"},
                ],
            }
        ]

        choice = inquire(questions)

        if choice == "A":
            upload_to_synology(args.project)

        elif choice == "D":
            ans = get_runs(args.project)
            if ans["status"]:
                run_paths = ans["runs"]
            else:
                continue
            
            #if not organize_trimmed(args.project, run_paths): continue
            if not organize_counts_by_asv(args.project): continue
            if not organize_counts_by_taxon(args.project): continue
            if not organize_logs(args.project, run_paths): continue
            if not organize_taxonomies_final(args.project): continue
            if not organize_taxonomies_intermediates(args.project): continue
            if not add_index(): continue
            if not add_references(args.project, run_paths): continue
            if not organize_reads(args.project, run_paths): continue
            if not share_unix_perms(args.project): continue

            if not remove_trash(args.project, run_paths): continue


        elif choice == "J":
            upload_to_jira(args.project)

        elif choice == "QUIT":
            sys.exit(0)
        elif choice == "B":
            args.project = choose_project()

    return


def ask_continue():
    if not confirm("Continue?"):
        return False
    return True

def ask_and_move(dirpath, filepaths):
    if len(filepaths) > 0:
        print()
        for filepath in filepaths:
            print(filepath)
        if confirm("Press y to move the above files:"):
            for filepath in filepaths:
                Path(filepath).rename(Path(dirpath) / Path(filepath))
        else:
            return False
    return True


# def archive(proj):
#     barcodes = []
#     for run in run_paths:
#         barcodes += glob.glob(str(run.resolve() / Path("barcodes.fastq")))
#     print(f"Found {len(barcodes)} barcode files in {len(run_paths)} runs.")

#     if(len(barcodes) > 0):
#         if not confirm("Delete?"):
#             return()

#         for barcode_filename in barcodes:
#             os.unlink(barcode_filename)

#     upload_to_synology(str(args.project))
#     return


# untested. takes one filepath string
def gz(filepath):
    with open(filepath, "rb") as f_in:
        with gzip.open(f"{filepath}.gz", "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
    Path(filepath).unlink()

    return


# untested. returns list of strings
def find_trimmed(run):
    files = glob.glob(str(run / Path("*.fastq")))
    return files


# untested. returns list of strings
def find_raw(run, ori):
    if ori == "fwd":
        return glob.glob(str(run / Path("fwdSplit/split_by_sample_out/*.fastq")))
    elif ori == "rev":
        return glob.glob(str(run / Path("revSplit/split_by_sample_out/*.fastq")))
    else:
        raise ValueError("Incorrect orientation. 'fwd' and 'rev' are allowed.")


# choose_upload_dest is a misnomer becaues it actually follows up by guiding user
# through the upload.
def upload_to_synology(dirpath):
    done = False
    while not done:
        questions = [
            {
                "type": "list",
                "name": "chooser",
                "message": "Upload project to Synology?",
                "choices": ["Upload", "Back"],
            }
        ]
        if inquire(questions) != "Upload":
            return

        user = getuser("Synology user: ")

        dirs_to_upload = [
            str(pathfinder['counts-by-ASV']),
            str(pathfinder['counts-by-taxon']),
            "MAP",
            "REPORT",
            "FASTQ",
            str(pathfinder['final_taxonomies']),
            str(pathfinder['intermediate_taxonomies']),
        ]

        try:
            synology.archive_project(
                user,
                dirpath.name,
                dirs_to_upload, 
            )
            done = True
        except Exception as e:
            print(e)

    return


def share_unix_perms(path):
    print("Expanding owner file permissions to group file permissions...")
    pathlist = Path(str(path)).rglob("*")
    user_to_group_divisor = int(stat.S_IRUSR / stat.S_IRGRP)
    for path in pathlist:
        try:
            st = path.stat()
        except Exception as e:
            print(str(e))
            if not ask_continue():
                return False
            else:
                continue
        usr_perms = st.st_mode & stat.S_IRWXU
        if usr_perms > 0:
            if os.geteuid() == st.st_uid:
                try:
                    path.chmod(st.st_mode | int(usr_perms / user_to_group_divisor))
                except PermissionError as e:
                    print(f"Unable to expand permissions to group for: {str(path)}\n{str(e)}")
                    # don't ask if the user wants to "continue". this is the last step, 
                    # so the user's choice doesn't make a difference
                    return False
            else:
                print(f"Cannot process {str(path)} without ownership\n")
                return False 
        else:
            print(f"Owner has no permissions on {str(path)}.\n")
            return False 

    return True


def getuser(prompt="User: "):
    try:
        default = os.getlogin()
    except OSError:
        default = os.environ["USER"]

    questions = [
        {"type": "input", "name": "user_prompt", "message": prompt, "default": default}
    ]
    return inquire(questions)


def upload_to_jira(proj):
    # CONNECT ##########################################################################
    options = {"server": "https://jira.igs.umaryland.edu/"}
    kerberos_options = {"mutual_authentication": "DISABLED"}

    connected = False
    while not connected:
        jira_user = getuser("Jira user:")
        jira_passwd = getpass.getpass(f"Password for {jira_user}:")
        basic_auth = (jira_user, jira_passwd)
        try:
            j_connection = JIRA(
                options=options,
                basic_auth=basic_auth,
                kerberos=True,
                kerberos_options=kerberos_options,
                max_retries=0,
            )
            connected = True
        except JIRAError as e:
            if e.status_code == 401:
                print("Unable to authenticate.")
                return
            else:
                print(e.response)
        except Exception as e:
            print(e)
            return

    # IDENTIFY TICKET ##################################################################
    issue_name = None
    while issue_name is None:
        questions = [
            {
                "type": "list",
                "name": "chooser",
                "message": "Find ticket to which to attach files:",
                "choices": [
                    {"name": "BROWSE OPEN MSL TICKETS", "value": "browse"},
                    {"name": "ENTER TICKET ID", "value": "enter"},
                    {"name": "BACK", "value": "back"},
                ],
            }
        ]
        choice = inquire(questions)

        if choice == "browse":
            print("Looking up issues...")
            issue_name = jira_choose_issue_from_list(j_connection)
        elif choice == "enter":
            issue_name = jira_find_issue_from_input(j_connection)
        else:
            return

    # UPLOAD ###########################################################################
    done = False
    while not done:
        files = sorted([str(child) for child in Path(".").iterdir() if child.is_file()])
        questions = [
            {
                "type": "checkbox",
                "name": "chooser",
                "message": "",
                "choices": [{"name": file} for file in files],
            }
        ]
        choices = inquire(questions, "Choose map(s) to upload:")
        if len(choices) == 0:
            if confirm("You selected 0 files. Continue?"):
                done = True
        else:
            for choice in choices:
                filepath = str(proj / Path(choice))
                if not attach_to_issue(j_connection, issue_name, filepath):
                    return
            done = True

    done = False
    while not done:
        files = sorted([str(child) for child in pathfinder["counts-by-taxon"].iterdir() if child.is_file()])
        # process = Popen(
        #     ["ls", "-al", str(proj / pathfinder["counts-by-taxon"])],
        #     stdin=PIPE,
        #     stdout=PIPE,
        #     stderr=PIPE,
        # )
        # tree = process.stdout.read().decode(encoding="utf8").split("\n")
        # if len(tree) < 5:
        #     print("Cannot read contents of current directory.")
        #     return
        questions = [
            {
                "type": "checkbox",
                "name": "chooser",
                "message": "",
                "choices": [{"name": file} for file in files],
            }
        ]
        choices = inquire(questions, "Choose count tables to upload:")

        if len(choices) == 0:
            if confirm("You selected 0 files. Continue?"):
                done = True
        else:
            for choice in choices:
                filepath = str(proj / pathfinder["counts-by-taxon"] / Path(choice))
                if not attach_to_issue(j_connection, issue_name, filepath):
                    return
            done = True

    return


def jira_choose_issue_from_list(j_connection):
    block_size = 100
    block_num = 0
    issues = {}
    while True:
        start_idx = block_num * block_size
        results = j_connection.search_issues(
            'project="MSL" AND (status="To Do" OR status="in progress") order by updated',
            start_idx,
            block_size,
        )
        if len(results) == 0:
            # Retrieve issues until there are no more to come
            break
        for res in results:
            issues[f"{res.key}: {res.fields.summary}"] = res.key
        block_num += 1

    questions = [
        {
            "type": "list",
            "name": "chooser",
            "message": "Choose ticket:",
            "choices": ["BACK"] + list(issues.keys()),
        }
    ]
    choice = inquire(questions)
    if choice != "BACK":
        return issues[choice]
    else:
        return None


def jira_find_issue_from_input(j_connection):
    # input prompt
    questions = [
        {
            "type": "input",
            "name": "user_prompt",
            "message": "Issue ID:",
            "default": "MSL-",
        }
    ]
    issue_name = inquire(questions)
    try:
        echo = j_connection.search_issues(f"id={issue_name}")
    except JIRAError as e:
        if e.status_code == 400:
            print("Unable to find issue.")
        else:
            print(e.text)
        return

    if len(echo) == 1:
        if echo[0].key == issue_name:
            if confirm(f"Use {echo[0].key}: {echo[0].fields.summary}?"):
                return issue_name
        else:
            print("Can't find issue.")
            return None
    else:
        print("Can't find issue.")
        return None


def attach_to_issue(j_connection, issue_name, filepath):
    file_path = Path(filepath).resolve()
    jira_issue = j_connection.issue(issue_name, expand="attachment")

    dupes = sum(
        [file_path.name == att.filename for att in jira_issue.fields.attachment]
    )

    if dupes > 0:
        prompt = (
            f"{dupes} file(s) named {file_path.name} already attached to {issue_name}."
        )

        questions = [
            {
                "type": "list",
                "name": "chooser",
                "message": prompt,
                "choices": [
                    "CANCEL",
                    "DELETE EXISTING/OVERWRITE",
                    "RENAME NEW ATTACHMENT",
                ],
            }
        ]

        choice = inquire(questions)

        if choice == "CANCEL":
            return False
        elif choice == "DELETE EXISTING/OVERWRITE":
            for att in jira_issue.fields.attachment:
                if file_path.name == att.filename:
                    j_connection.delete_attachment(att.id)

        elif choice == "RENAME NEW ATTACHMENT":
            questions = [
                {
                    "type": "input",
                    "name": "user_prompt",
                    "message": "New filename:",
                    "default": f"{file_path.name}.1",
                }
            ]
            with tempfile.TemporaryDirectory() as tmpdirname:
                new_file_path = Path(tmpdirname) / Path(inquire(questions))
                shutil.copy(file_path, new_file_path)
                with new_file_path.open(mode="rb") as f:
                    j_connection.add_attachment(
                        issue_name, f
                    )  # Please check if attachment already exists, and if so, ask to delete
                print(
                    f"Uploaded {str(file_path.relative_to(args.project))} to {issue_name} as {new_file_path.name}"
                )
                return True

    with file_path.open(mode="rb") as f:
        j_connection.add_attachment(
            issue_name, f
        )  # Please check if attachment already exists, and if so, ask to delete
    print(f"Uploaded {str(file_path.relative_to(args.project))} to {issue_name}")
    return True


def inquire(questions, extra_prompt=None):
    print("\n")
    if extra_prompt is not None:
        print(extra_prompt)
    answers = PyInquirer.prompt(questions, style=examples.custom_style_2)
    if answers:
        for answer in answers.values():
            return answer  # only one answer expected
    else:
        sys.exit(0)


def choose_project(default=None):
    def helper():
        # simple browser, chooses from current working directory with no navigation allowed
        if confirm("Is the project in the current directory?"):
            return Path.cwd()
        else:
            while True:
                process = Popen(
                    ["ls", "-al", str(Path.cwd())], stdin=PIPE, stdout=PIPE, stderr=PIPE
                )
                tree = process.stdout.read().decode(encoding="utf8").split("\n")
                if len(tree) < 5:
                    print("Cannot read contents of current directory. Exiting...")
                    sys.exit(0)
                questions = [
                    {
                        "type": "list",
                        "name": "chooser",
                        "message": "Please choose project directory.",
                        "choices": ["QUIT"] + tree[3 : len(tree) - 1],
                    }
                ]
                choice = inquire(questions)

                if choice[0] == "d":
                    basename = choice.split()[8]
                    fullname = (Path.cwd() / Path(basename)).resolve()
                    return fullname
                elif choice[0] == "l":
                    target = Path(choice.split()[8]).resolve()
                    if target.is_dir():
                        return target
                elif choice == "QUIT":
                    sys.exit(0)

                print("\nNot a directory!\n")

    if default is None:
        project = helper()
    else:
        try:
            project = Path(default).resolve()  # May raise FileNotFoundError
        except FileNotFoundError as e:
            print(e)
            sys.exit(1)
    os.chdir(str(project))

    return project


def confirm(prompt, default=False):
    questions = [
        {"type": "confirm", "name": "response", "message": prompt, "default": default}
    ]
    return inquire(questions)

def get_runs(proj_path):
    ans = {
        "status": True,
        "runs": []
    }
    subdirs = [x for x in Path(proj_path).iterdir() if x.is_dir()]

    if len(subdirs) == 0:
        if confirm("Project contains no subdirectories. Continue?"):
            return []
        else:
            return ()
    choices = []

    for subdir in subdirs:
        choice = {"name": subdir.name}
        if subdir.name in list(project_metadata()['runs'].keys()):
            choice["checked"] = True
        else:
            choice["checked"] = False
        choices.append(choice)
    choices.sort(key=(lambda choice: choice["name"]))

    prompt = "Use <SPACE> to choose runs. Press <ENTER> to confirm."
    questions = [
        {"type": "checkbox", "name": "chooser", "message": "", "choices": choices}
    ]

    response = inquire(questions, prompt)
    if len(response) == 0:
        if confirm("You selected 0 runs. Continue?"):
            pass
        else:
            ans["status"] = False
    else:
        ans["runs"] = [proj_path / Path(subdir) for subdir in response]

    return ans

def project_metadata():
    try:
        with Path(".meta.json").open("r") as fh:
            metadata = json.load(fh)
    except FileNotFoundError:
        raise Exception("Cannot find project metadata from Part 2. Re-run Part 2.")
    if not type(metadata) is dict:
        raise Exception("Project metadata not a JSON dictionary")

    return metadata

def read_metadata(run_path):
    with (run_path / Path(".meta.json")).open("r") as fh:
        run_info = json.load(fh)
    if not type(run_info) is dict:
        raise Exception("Run metadata not a JSON dictionary")
    if not "params" in run_info.keys() or not "checkpoints" in run_info.keys():
        raise Exception("Run metadata doesn't contain 'checkpoints' and 'params'")
    if not type(run_info["params"]) is dict or not type(run_info["checkpoints"]) is dict:
        raise Exception("Run metadata 'params' or 'checkpoints' is not a JSON dictionary")

    return run_info

def organize_reads(proj_path, run_paths):
    try:
        organized_dir = "FASTQ"
        filepaths = []
        any_files = False

        executor = CONFIG['executor']
        external_exec = False
        print("\nThis pipeline's default executor is:\n\n" + executor)
        if confirm("Use executor to transfer large files? (Choose yes if not using screen or tmux)", default=True):
            external_exec = True

        for run_path in run_paths:
            ill_fwd_dir = proj_path / run_path / Path("fwdSplit") / Path("split_by_sample_out")
            ill_rev_dir = proj_path / run_path / Path("revSplit") / Path("split_by_sample_out")

            # try:
            #     run_info = read_metadata(proj_path / run_path)
            # except Exception as e:
            #     print(str(e))
            #     print("WARNING: Can't open run metadata from " + run_path.name + ". Unable to locate raw reads")
            #     return True # NOT A SHOWSTOPPER

            # try:
            #     filepaths += [
            #         os.path.join(run_path, rel_path) for rel_path in run_info['checkpoints']["samples"].keys()
            #     ]
            # except KeyError:
            #     print("WARNING: Incompatible with the pipeline version used on this run. Cannot find raw read files.")
            #     return True # not a show-stopper

            if ill_fwd_dir.is_dir() and ill_rev_dir.is_dir():
                nFiles = len(list(ill_fwd_dir.iterdir()))
                nFiles += len(list(ill_rev_dir.iterdir()))
                any_files = True
                
                make_for_contents(organized_dir, any_files)
                subdir_destination = str(Path(organized_dir) / Path(run_path.name))
                make_for_contents(subdir_destination)

                if external_exec:
                    print(f"Using executor to copy {nFiles} raw read files to {subdir_destination}.\n")
                    
                    fwd_wildcard = str(ill_fwd_dir / Path("*"))
                    rev_wildcard = str(ill_rev_dir / Path("*"))

                    cmd = "module load sge && " + executor + " -V \'cp -f %s " + subdir_destination + "\'"
                    print(cmd % fwd_wildcard)
                    run(cmd % fwd_wildcard, shell=True)
                    print(cmd % rev_wildcard)
                    run(cmd % rev_wildcard, shell=True)

                    print("\nRemember to delete any stray STDOUT and STDERR files from the working directory.\n")

                    time.sleep(5)

                else:
                    print(f"Copying {nFiles} raw read files to {subdir_destination}.")
                    for filepath in filepaths:
                        shutil.copy2(filepath, Path(subdir_destination) / Path(filepath).name)

                # gzip all if needed
                for f in Path(subdir_destination).iterdir():
                    if f.suffix != ".gz":
                        gz(str(f))
                            

        if not any_files:
            print(f"Skipping creation of {organized_dir} (no demuxed reads found)")
    except Exception as e:
        print(str(e))
        if not ask_continue():
            return False

    return True


# def organize_trimmed(proj_path, run_paths):
#     try:
#         organized_dir = "FASTQ_TRIMMED"

#         filepaths = []
#         any_files = False
#         for run_path in run_paths:
#             filepaths += glob.glob(str(proj_path / run_path / Path("*_tc.fastq*")))
#             filepaths += glob.glob(
#                 str(proj_path / run_path / Path("tagcleaned") / Path("*.fastq*"))
#             )

#             if len(filepaths) > 0:
#                 make_for_contents(organized_dir, any_files)
#                 any_files = True
#                 subdir_destination = str(Path(organized_dir) / run_path.name)
#                 make_for_contents(subdir_destination)
#                 print(
#                     f"Copying {len(filepaths)} adapter-trimmed read files to {subdir_destination}."
#                 )
#                 for filepath in filepaths:
#                     shutil.copy2(filepath, Path(subdir_destination) / Path(filepath).name)

#                 for f in Path(subdir_destination).iterdir():  # gzip if necessary
#                     if f.suffix != ".gz":
#                         gz(str(f))

#         if not any_files:
#             print(f"Skipping creation of {organized_dir} (no adapter-trimmed reads found)")

#     except Exception as e:
#         print(str(e))
#         if not ask_continue():
#             return False
    
#     return True


def organize_counts_by_asv(proj_path):
    try:
        organized_dir = str(pathfinder['counts-by-ASV'])

        filepaths = glob.glob(str(proj_path / Path("*all_runs_dada2_abundance_table.rds")))
        filepaths += glob.glob(str(proj_path / Path("*all_runs_dada2_abundance_table.csv")))
        filepaths += glob.glob(str(proj_path / Path("*asvs+taxa.csv")))

        if len(filepaths) > 0:
            make_for_contents(organized_dir)
            print(f"Moving {len(filepaths)} sample-ASV count tables to {organized_dir}/.")
            for filepath in filepaths:
                Path(filepath).rename(Path(organized_dir) / Path(filepath).name)
        else:
            print(f"Skipping creation of {organized_dir} (no count-by-ASV tables found)")

    except Exception as e:
        print(str(e))
        if not ask_continue():
            return False

    return True


def organize_counts_by_taxon(proj_path):
    try:
        organized_dir = str(pathfinder['counts-by-taxon'])

        filepaths = glob.glob(str(proj_path / Path("*taxa-merged.csv")))

        if len(filepaths) > 0:
            make_for_contents(organized_dir)
            print(f"Moving {len(filepaths)} sample-taxon count tables to {organized_dir}/.")
            for filepath in filepaths:
                Path(filepath).rename(Path(organized_dir) / Path(filepath).name)
        else:
            print(f"Skipping creation of {organized_dir} (no count-by-taxon tables found)")

    except Exception as e:
        print(str(e))
        if not ask_continue():
            return False

    return True


def make_for_contents(dir, silent=False):
    try:
        os.makedirs(dir)
        print(f"Creating ./{dir}/")
    except FileExistsError:
        if not silent:
            print(f"Directory ./{dir}/ already exists")
    except Exception as e:
        raise e

def organize_logs(proj_path, run_paths):
    try:
        organized_dir = "LOGS"

        # Part 1 logs found in the run directories. Just copy.

        filepaths = glob.glob(str(proj_path / Path("*pipeline_log.txt")))
        for run in run_paths:
            filepaths += glob.glob(str(proj_path / run / Path("*pipeline_log.txt")))

        if len(filepaths) > 0:
            make_for_contents(organized_dir)
            print(f"Moving {len(filepaths)} log files to LOGS/.")

            for filepath in filepaths:
                shutil.copy(Path(filepath), Path(organized_dir) / Path(filepath).name)
        else:
            print(f"Skipping creation of {organized_dir} (no logs found)")

    except Exception as e:
        print(str(e))
        if not ask_continue():
            return False
    
    return True


def organize_taxonomies_final(proj_path):
    try:
        organized_dir = str(pathfinder['final_taxonomies'])

        filepaths = glob.glob(str(proj_path / "cmb_tx.txt"))

        if len(filepaths) > 0:
            make_for_contents(organized_dir)
            print(f"Moving cmb_tx.txt to {organized_dir}/.")

            for filepath in filepaths:
                Path(filepath).rename(Path(organized_dir) / Path("final_taxonomy.txt"))
        else:
            print(f"Skipping creation of {organized_dir} (no final taxonomy files found)")

    except Exception as e:
        print(str(e))
        if not ask_continue():
            return False
    return True


def organize_taxonomies_intermediates(proj_path):
    try:
        organized_dir = str(pathfinder['intermediate_taxonomies'])
        any_files = False

        filepaths = glob.glob(str(proj_path / Path("MC_order7_results.txt")))
        if len(filepaths) > 0:
            any_files = True
            make_for_contents(organized_dir)
            print(f"Moving MC_order7_results.txt to {organized_dir}/PECAN_raw.txt.")

            for filepath in filepaths:
                Path(filepath).rename(Path(organized_dir) / Path("PECAN_raw.txt"))

        filepaths = glob.glob(str(proj_path / Path("*SILVA*.classification.csv")))
        filepaths += glob.glob(str(proj_path / Path("*HOMD*.classification.csv")))
        filepaths += glob.glob(str(proj_path / Path("*UNITE*.classification.csv")))

        if len(filepaths) > 0:
            any_files = True
            make_for_contents(organized_dir, any_files)
            print(f"Moving {len(filepaths)} raw RDP classifier files to {organized_dir}/")

            for filepath in filepaths:
                newname = re.sub(
                    r"^.*?([^_]*)\.classification\.csv$", r"\1_raw.csv", Path(filepath).name
                )
                Path(filepath).rename(
                    Path(organized_dir) / Path(newname)
                )

        filepaths = glob.glob(str(proj_path / Path("silva_condensed.txt")))
        filepaths += glob.glob(str(proj_path / Path("homd_condensed.txt")))
        filepaths += glob.glob(str(proj_path / Path("unite_condensed.txt")))
        if len(filepaths) > 0:
            make_for_contents(organized_dir, any_files)
            print(f"Moving {len(filepaths)} condensed taxonomy file to {organized_dir}/.")

            for filepath in filepaths:
                match = re.match(r"^([^_]*)_condensed\.txt$", Path(filepath).name)
                tax_name = match.group(1)
                Path(filepath).rename(Path(organized_dir) / Path(f"{tax_name.upper()}.txt"))

        if not any_files:
            print(
                f"Skipping creation of {organized_dir} (no intermediate taxonomy files found)"
            )

    except Exception as e:
        print(str(e))
        if not ask_continue():
            return False
    
    return True

def add_index():
    try:
        with open("_INDEX.txt", "w") as outfile:
            outfile.write(index_contents)
    except Exception as e:
        print(str(e))
        if not ask_continue():
            return False
    return True

def add_references(proj_path, run_paths):
    # get illumina references as orderedDict
    references = Citations()
    print('')

    try: # ugh I don't like this huge try catch block, but it works...and I don't
        # know what exceptions to catch any way
        # add starter references for illumina or pacbio runs
        run_types = set()
        for run_path in run_paths:
            try:
                run_info = read_metadata(proj_path / run_path)
                run_type = run_info['params']['platform'].upper()
                if run_type not in CITATION_GROUPS.keys():
                    print(f"Unrecognized run type: {run_type}. Don't know what citations to add.")
                else:
                    run_types.add(run_type)
            except Exception as e:
                print(f"WARNING: Unable to determine what platform was used in {run_path.name}. Cannot add required citations.")
            
        if len(run_types) > 0:
            print(f"Adding citations for the following run types:")
            for run_type in run_types:
                print("\t" + run_type)
                references.merge(CITATION_GROUPS[run_type])
            print("")
        else:
            print(f"No recognized run types found.\n")
        
        # add citations for any taxonomy/classifiers used
        directory = str(pathfinder['intermediate_taxonomies'])
        taxonomy_files = glob.glob(os.path.join(directory, "*_raw.*"))
        for my_file in taxonomy_files:
            match = re.match('^(.*?)_raw.(?:csv|txt)$', Path(my_file).name)
            if match:
                taxonomy = match.group(1)
                if taxonomy not in CITATION_GROUPS.keys():
                    print(f"Unrecognized taxonomy: {taxonomy} inferred from file {my_file}, unable to add citations.")
                    continue
                references.merge(CITATION_GROUPS[taxonomy])
                print(f"Adding citation for {taxonomy}")

        # write to file
        with open("_REFERENCES.txt", "w") as outfile:
            outfile.write(references.to_string())

    except Exception as e:
        print(e)
        if not ask_continue():
            return False
    return True

def remove_trash(proj_path, run_paths):
    try:
        trash_files = [] # if we need to use os.unlink
        trash_dirs = [] # if we need to use shutil.rmtree
        for run_path in run_paths:
            if run_path.is_symlink():
                for fastq in run_path.glob('*Split/seqs.fastq'):
                    trash_files.append(str(fastq))
                for fastq in run_path.glob('barcodes.fastq'):
                    trash_files.append(str(fastq))
                for fastq in run_path.glob('*_tc.fastq.gz'):
                    trash_files.append(str(fastq))
                
                # we must delete the symlink or directory last
                trash_files.append(
                    {
                        'name': str(run_path) + " (L)", 
                        'value': str(run_path)
                    }
                )
            else:
                trash_dirs.append(str(run_path))

        for entry in os.scandir("."):
            if (
                entry.name.startswith(".")
                or entry.name.startswith("rTmp.")
                or entry.name.startswith("dump.rda")
                or re.match(".*\.[eo][0-9]+$", entry.name)
            ) and entry.is_file():
                trash_files.append(entry.name)

        if len(trash_files) > 0 or len(trash_dirs) > 0:
            print()
            print("\n".join(trash_dirs))
            trash_files_names = [x["name"] if isinstance(x, dict) else x for x in trash_files]
            print("\n".join(trash_files_names))

            questions = [
                {
                    "type": "list",
                    "name": "chooser",
                    "message": "The above files/directories will be deleted. (L) indicates removal of a symlink, while the real files will be preserved. Ensure that all important result files have been copied out.",
                    "choices": ["Keep", "Delete all"],
                }
            ]
            if not inquire(questions) == "Delete all":
                return True
            
            # part2_danger = not (pathfinder['counts-by-taxon'].is_dir() and pathfinder['counts-by-ASV'].is_dir() and pathfinder['intermediate_taxonomies'].is_dir() and pathfinder['final_taxonomies'].is_dir())
            confirm_msg = "Are you sure?"
            # if part2_danger:
            #     confirm_msg = "Looks like Part 2 hasn't been run!  " + confirm_msg
            questions = [
                {
                    "type": "list",
                    "name": "chooser",
                    "message": confirm_msg,
                    "choices": ["Keep", "Delete all"],
                }
            ]
            if not inquire(questions) == "Delete all":
                return True
            for trash_dir in trash_dirs:
                shutil.rmtree(trash_dir)
            for trash_file in trash_files:
                trash_file = trash_file['value'] if isinstance(trash_file, dict) else trash_file
                os.unlink(trash_file)

        return True
    except Exception as e:
        print(str(e))
        if not ask_continue():
            return False

    return True

# def share_project(proj_path):
#     try:
#         run(["sharedir", "./"])

#         print("Read and write permissions granted for group.")
#         return
#     except Exception as e:
#         print(str(e))
#         return False

index_contents = """
*_all_runs_dada2_ASV.fasta: Denoised ASVs (Amplicon Sequence Variants).
*_DADA2_stats.txt: Table of sample read throughput.
COUNTS_BY_ASV: Read count tables with reads counted by ASV
COUNTS_BY_TAXON: Read count tables with reads counted by taxon. Corresponding taxonomic assignment files in ./TAXONOMY_FINAL and ./TAXONOMY_INTERMEDIATES.
FASTQ: Raw reads per sample
FASTQ_TRIMMED: Sample reads trimmed of primers
LOGS: Documentation of MSL computational pipeline parameters and workflow.
MAPS: Barcodes used to demultiplex reads to samples. This directory may be absent if all runs were PacBio runs, or if files sent to MSL did not require demultiplexing.
REPORT: An HTML document and accompanying assets summarizing the results. To display the report in a browser, all contents of the directory must be present.
TAXONOMY_FINAL: File(s) documenting taxonomies assigned to ASVs, drawing on results in TAXONOMY_INTERMEDIATES. Each taxonomically-annotated file in ./COUNTS_BY_ASV and ./COUNTS_BY_TAXON is named after a file in this directory.
TAXONOMY_INTERMEDIATES: Taxonomic classifications rendered by one or more taxonomic classifiers, including full taxonomic output from the RDP classifier.
"""

class Step_citations(set):
    def __init__(self, title, citations=[], order=None):
        if order is not None:
            try:
                self.order = int(order)
            except:
                raise TypeError("Cannot cast Step_citations order to int")
        else:
            self.order = order
        try:
            self.title = str(title)
        except:
            raise TypeError("Cannot cast Step_citations title to str")
        try:
            self.update([str(citation) for citation in citations])
        except:
            raise TypeError("One or more citations cannot be cast to str")
    
    def to_string(self):
        return("\n".join(self))
    
    def merge(self, other):
        if type(other) != Step_citations:
            raise TypeError("Attempt to merge a Step_citations object with an incompatible object.")
        if self.title != other.title:
            raise Exception(f"Not allowed to merge citations for different steps: {self.title} and {other.title}")
        if self.order and other.order and self.order != other.order:
            raise Exception(f"Conflicting order of citation step {self.title}: {self.order} or {other.order}?")
        if (not self.order) and other.order:
            self.order = other.order
        self.update(other)


class Citations(dict):
    def __init__(self, steps=[]):
        for step in steps:
            if type(step) is not Step_citations:
                raise TypeError("Citations can only hold Step_citations objects")
            self[step.title] = step
    def merge(self, other):
        if type(other) != Citations:
            raise TypeError("Attempt to merge a Citations object with an incompatible object.")
        for step_title, step in other.items():
            if step_title in self.keys():
                self[step_title].merge(other[step_title])
            else:
                self[step_title] = step
    
    def to_string(self):
        sorted_steps = sorted(self.keys(), key=lambda step_title: self[step_title].order)
        content = ""
        for step in sorted_steps:
            content += self[step].title + ":\n"
            content += self[step].to_string() + "\n\n"
        return(content)

CITATION_GROUPS = {
    'ILLUMINA': Citations([
        Step_citations(
            title='Barcode extraction, library demultiplexing', 
            order=1,
            citations=[
                'QIIME allows analysis of high-throughput community sequencing data. J Gregory Caporaso, Justin Kuczynski, Jesse Stombaugh, Kyle Bittinger, Frederic D Bushman, Elizabeth K Costello, Noah Fierer, Antonio Gonzalez Pena, Julia K Goodrich, Jeffrey I Gordon, Gavin A Huttley, Scott T Kelley, Dan Knights, Jeremy E Koenig, Ruth E Ley, Catherine A Lozupone, Daniel McDonald, Brian D Muegge, Meg Pirrung, Jens Reeder, Joel R Sevinsky, Peter J Turnbaugh, William A Walters, Jeremy Widmann, Tanya Yatsunenko, Jesse Zaneveld and Rob Knight; Nature Methods, 2010; doi:10.1038/nmeth.f.303'
            ]
        ),
        Step_citations(
            'Primer removal',
            [
                'Schmieder R, Lim YW, Rohwer F, Edwards R: TagCleaner: Identification and removal of tag sequences from genomic and metagenomic datasets. BMC Bioinformatics 2010, 11:341. [PMID: 20573248]'
            ],
            2
        ),
        Step_citations(
            'QC, decontamination, and denoising',
            [
                'Callahan, B., McMurdie, P., Rosen, M. et al. DADA2: High-resolution sample inference from Illumina amplicon data. Nat Methods 13, 581–583 (2016). https://doi.org/10.1038/nmeth.3869'
            ],
            3
        ),
        Step_citations(
            'Taxonomic assignment',
            [
                'Callahan, B., McMurdie, P., Rosen, M. et al. DADA2: High-resolution sample inference from Illumina amplicon data. Nat Methods 13, 581–583 (2016). https://doi.org/10.1038/nmeth.3869',
                "Wang, Q, G. M. Garrity, J. M. Tiedje, and J. R. Cole. 2007. Naïve Bayesian Classifier for Rapid Assignment of rRNA Sequences into the New Bacterial Taxonomy. Appl Environ Microbiol. 73(16):5261-7."
            ],
            10
        )
    ]),
    'PACBIO': Citations([
        Step_citations(
            'QC, decontamination, and denoising',
            [
                'Callahan, B., McMurdie, P., Rosen, M. et al. DADA2: High-resolution sample inference from Illumina amplicon data. Nat Methods 13, 581–583 (2016). https://doi.org/10.1038/nmeth.3869'
            ],
            3
        ),
        Step_citations(
            'Taxonomic assignment',
            [
                "Wang, Q, G. M. Garrity, J. M. Tiedje, and J. R. Cole. 2007. Naïve Bayesian Classifier for Rapid Assignment of rRNA Sequences into the New Bacterial Taxonomy. Appl Environ Microbiol. 73(16):5261-7."
            ],
            10
        ),
    ]),
    'SILVA132': Citations([
        Step_citations(
            'Taxonomic assignment',
            [
                'Quast C, Pruesse E, Yilmaz P, Gerken J, Schweer T, Yarza P, Peplies J, Glöckner FO (2013) The SILVA ribosomal RNA gene database project: improved data processing and web-based tools. Nucl. Acids Res. 41 (D1): D590-D596.',
                'Yilmaz P, Parfrey LW, Yarza P, Gerken J, Pruesse E, Quast C, Schweer T, Peplies J, Ludwig W, Glöckner FO (2014) The SILVA and "All-species Living Tree Project (LTP)" taxonomic frameworks. Nucl. Acids Res. 42:D643-D648'
                'Callahan, Benjamin. (2018). Silva taxonomic training data formatted for DADA2 (Silva version 132) [Data set]. Zenodo.'
            ]
        )
    ]),
    'SILVA138': Citations([
        Step_citations(
            'Taxonomic assignment',
            [
                'Quast C, Pruesse E, Yilmaz P, Gerken J, Schweer T, Yarza P, Peplies J, Glöckner FO (2013) The SILVA ribosomal RNA gene database project: improved data processing and web-based tools. Nucl. Acids Res. 41 (D1): D590-D596.',
                'Yilmaz P, Parfrey LW, Yarza P, Gerken J, Pruesse E, Quast C, Schweer T, Peplies J, Ludwig W, Glöckner FO (2014) The SILVA and "All-species Living Tree Project (LTP)" taxonomic frameworks. Nucl. Acids Res. 42:D643-D648'
                'Michael R. McLaren. (2020). Silva SSU taxonomic training data formatted for DADA2 (Silva version 138) (Version 1) [Data set]. Zenodo. http://doi.org/10.5281/zenodo.3986799'
            ]
        )
    ]),
    'SILVA138forPB': Citations([
        Step_citations(
            'Taxonomic assignment',
            [
                'Quast C, Pruesse E, Yilmaz P, Gerken J, Schweer T, Yarza P, Peplies J, Glöckner FO (2013) The SILVA ribosomal RNA gene database project: improved data processing and web-based tools. Nucl. Acids Res. 41 (D1): D590-D596.',
                'Yilmaz P, Parfrey LW, Yarza P, Gerken J, Pruesse E, Quast C, Schweer T, Peplies J, Ludwig W, Glöckner FO (2014) The SILVA and "All-species Living Tree Project (LTP)" taxonomic frameworks. Nucl. Acids Res. 42:D643-D648'
            ]
        )
    ]),
    'PECAN': Citations([
        Step_citations(
            'Taxonomic assignment',
            [
                'Pawel Gajer, Jacques Ravel, Johanna Holm. Github. 2018. SpeciateIT. [Online]. Available: https://github.com/Ravel-Laboratory/speciateIT'
            ]
        )
    ]),
    'UNITE': Citations([
        Step_citations(
            'Taxonomic assignment',
            [
                'UNITE Community (2017): UNITE general FASTA release. Version 01.12.2017. UNITE Community. https://doi.org/10.15156/BIO/587475'
            ]
        )
    ]),
    'HOMD': Citations([
        Step_citations(
            'Taxonomic assignment',
            [
                'F. Escapa, I., Huang, Y., Chen, T., Lin, M., Kokaras, A., Dewhirst F.E., Lemon, K.P. (2020) Construction of habitat-specific training sets to achieve species-level assignment in 16S rRNA gene datasets. Microbiome 8, 65. Online Open Access https://doi.org/10.1186/s40168-020-00841-w'
            ]
        )
    ])
}

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
    touch authorized_keys
    cat ~/id_rsa.pub >> ~/.ssh/authorized_keys
    rm ~/id_rsa.pub
    chmod 600 authorized_keys

    Ensure this command works:
    rsync --list-only SYNOLOGY_USER@10.90.233.174:/volume1/NetBackup/MSL


            """
    )

    parser.add_argument("project", nargs="?", help="Project folder", default=None)

    args = parser.parse_args()

    main(args)
