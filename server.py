import os
import shutil
from flask import Flask, request
from flask_cors import CORS

import default_path

print("Hello World")

uploads_reads_align_const = default_path.upload_reads_align
os.makedirs(uploads_reads_align_const, exist_ok=True)
uploads_refs_align_const = default_path.upload_refs_align
os.makedirs(uploads_refs_align_const, exist_ok=True)
output_dir_align = default_path.output_dir_align
print(output_dir_align)
os.makedirs(output_dir_align, exist_ok=True)

# Setup flask server
app = Flask(__name__)
CORS(app)


@app.route('/align', methods=['GET', 'POST'])
# Route to run primalscheme multiplexer on a fasta file
def align():

    if len(os.listdir(uploads_reads_align_const)) != 0:
        EraseFile(uploads_reads_align_const)
        print("Good1")
    if len(os.listdir(uploads_refs_align_const)) != 0:
        EraseFile(uploads_refs_align_const)
        print("Good2")
    # if len(os.path(output_dir_align)) != 0:
    #     shutil.rmtree(output_dir_align)
    #     print("Good3")

    run_name = request.form["NameOutput"]
    output_dir_run = f"{output_dir_align}{run_name}"
    os.makedirs(output_dir_run, exist_ok=True)
    print(run_name)

    # Load fastq zip
    fastq = request.files["Read"]
    fastq.save(os.path.join(uploads_reads_align_const, fastq.filename))
    print("FastQ saved")

    # Load list of files name for references
    tmp_files = request.form["Refs"]
    list_files = list(tmp_files.split(","))
    # Save each file
    for i in range(len(list_files)):
        file = request.files[list_files[i]]
        print(file.filename)
        file.save(os.path.join(uploads_refs_align_const, file.filename))
        print(f"Fasta file {file.filename} saved")

    print("trop chelou!!")
    # os.system("nextflow run execute.nf")

    return "Done"


def EraseFile(path):
    files = os.listdir(path)
    for i in range(0, len(files)):
        os.remove(path + files[i])


if __name__ == "__main__":
    app.run(port=8080)
