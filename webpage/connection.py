# Here's an example program in Python that will accept a file and connect to a remote server via ssh key, execute a shell script,
# and download the resulting output files to the local computer using the 'paramiko' SSH library:
"""
Steps:

0. prepare config and input files (fasta QC here again sprotein or nuclaotide)
1. send to remote server input params and input file
2. execute pipeline with params on the server in custom dir
3. zippify output and prepare image for interface
4. send files to local
DONE 5. send email to user with zip file
6. store data in the local database (just some results dir and subdirs that contain results for each output)
7. access the files from results page in the interface (easy read files)
8. remove files after X hours (24 or 48) (run some process that will check dir creation time and delete if its life excess X h)

"""


import random
import os
import tempfile
from typing import List

import email, smtplib, ssl

from email import encoders
from email.mime.base import MIMEBase
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText

import paramiko

def connect():
    # Define ssh connection parameters
    ssh_host = os.environ['REMOTE_HOST']
    ssh_port = os.environ['REMOTE_PORT']
    ssh_user = os.environ['REMOTE_USER']
    ssh_key  = os.environ['REMOTE_KEY']

    # Define local directory to download files to
    local_dir = '.'

    # Define remote directory where the script and output files are located
    remote_dir = tempfile.mktemp('', 'results_', 'nemu_workdir')

    # Define shell script to execute on remote server
    script = 'mkdir output_dir && screenfetch -Nn > output_dir/some_extra_file.txt'
    script = 'nextflow run -c ...'

    # Create ssh client and connect to remote server
    ssh = paramiko.SSHClient()
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    ssh.connect(ssh_host, ssh_port, username=ssh_user, key_filename=ssh_key)

    # Execute script on remote server
    stdin, stdout, stderr = ssh.exec_command(script)
    exit_status = stdout.channel.recv_exit_status()

    if exit_status != 0:  # print error message if script exits with non-zero status
        print(stderr.readlines())
    else:  # download output files if script executes successfully
        sftp = ssh.open_sftp()
        remote_dir_full = os.path.join(remote_dir, 'output_dir')
        remote_files = sftp.listdir(remote_dir_full)  # list output files in remote directory
        for remote_file in remote_files:
            sftp.get('%s/output_dir/%s' % (remote_dir, remote_file), os.path.join(local_dir, remote_file))  # download each output file to local_dir
        sftp.close()


    # Close ssh client
    ssh.close()


def send_email(receiver_email: str, filenames: List[str], subject=None, body=None):
    """
    params:
        - filename: str - path to filename in same directory as script
    """
    subject = subject or "NeMu pipeline execution results"
    body = body or "This is an email with attachment sent from Python"
    sender_email = os.environ['EMAIL_LOGIN']
    password = os.environ['EMAIL_PASSWORD']

    # Create a multipart message and set headers
    message = MIMEMultipart()
    message["From"] = sender_email
    message["To"] = receiver_email
    message["Subject"] = subject
    message["Bcc"] = receiver_email  # Recommended for mass emails

    # Add body to email
    message.attach(MIMEText(body, "plain"))

    # Open PDF file in binary mode
    for filename in filenames:
        with open(filename, "rb") as attachment:
            # Add file as application/octet-stream
            # Email client can usually download this automatically as attachment
            part = MIMEBase("application", "octet-stream")
            part.set_payload(attachment.read())

        # Encode file in ASCII characters to send by email    
        encoders.encode_base64(part)
        # Add header as key/value pair to attachment part
        part.add_header(
            "Content-Disposition",
            f"attachment; filename= {os.path.basename(filename)}",
        )
        # Add attachment to message and convert message to string
        message.attach(part)

    text = message.as_string()

    # Log in to server using secure context and send email
    context = ssl.create_default_context()
    with smtplib.SMTP_SSL("smtp.mail.ru", 465, context=context) as server:
        server.login(sender_email, password)
        server.sendmail(sender_email, receiver_email, text)


if __name__ == "__main__":
    send_email('galogenid36326@gmail.com', ['requirements.txt', '.streamlit/config.toml'])



# To use this program, save it to a file (e.g. remote_execution.py) and modify the ssh_host, 
# ssh_username, ssh_key_file, local_dir, remote_dir, and script variables to match your specific use case. 
# Then, run the program from the command line using:

# python remote_execution.py /path/to/file


# where /path/to/file is the path to the file you want to use as an input argument. 
# The program will execute the specified shell script on the remote server, and download 
# any resulting output files to the specified local directory