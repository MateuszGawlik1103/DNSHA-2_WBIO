from flask import Flask, request, render_template
from dnsha512 import sha512_hash, add_padding, text_to_nuc_list

app = Flask(__name__)

@app.route('/', methods=['GET', 'POST'])
def home():
    hash_result = None
    input_text = None
    if request.method == 'POST':
        input_text = request.form['input_text']
        nuc_list = [add_padding(text_to_nuc_list(input_text))]
        dna_hash = sha512_hash(nuc_list)

        hash_result = str(dna_hash)
        hash_result = hash_result.replace('], ', '],<br>')
    return render_template('index.html', input_text=input_text, hash_result=hash_result)

if __name__ == '__main__':
    app.run(debug=True)
