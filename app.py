from flask import Flask, request, render_template, redirect, url_for
from dnsha512 import dnsha512_hash, hash_photo, add_padding, text_to_nuc_list, binary_to_nuc_list
import os

app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = 'uploads'
app.config['MAX_CONTENT_LENGTH'] = 16 * 1024 * 1024  # 16 MB limit

# Ensure the upload folder exists
os.makedirs(app.config['UPLOAD_FOLDER'], exist_ok=True)

@app.route('/', methods=['GET', 'POST'])
def home():
    hash_result = None
    input_text = None
    if request.method == 'POST':
        if 'input_text' in request.form:
            input_text = request.form['input_text']
            nuc_list = [add_padding(text_to_nuc_list(input_text))]
            dna_hash = dnsha512_hash(nuc_list)

            hash_result = str(dna_hash)
            hash_result = hash_result.replace('], ', '],<br>')

        if 'file' in request.files:
            file = request.files['file']
            if file and file.filename != '':
                file_path = os.path.join(app.config['UPLOAD_FOLDER'], file.filename)
                file.save(file_path)
                hash_result = hash_photo(file_path)
                hash_result = str(hash_result)
                hash_result = hash_result.replace('], ', '],<br>')

    return render_template('index.html', input_text=input_text, hash_result=hash_result)

if __name__ == '__main__':
    app.run(debug=True)
