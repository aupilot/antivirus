# from os import mkdir, path
import numpy as np
from transformers import T5EncoderModel, T5Tokenizer, BertModel, BertTokenizer
import torch
from Bio.Seq import Seq
import sentencepiece            # required!
from scipy.spatial import distance
from scipy.fftpack import dct, idct
from helper import highlight_differences

alphabet = [
    "_", "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y",
]

# pip install transformers sentencepiece

class EmbedT5(object):
    def __init__(self, reference):
        # weights
        # if not path.exists("prot_t5"):
        #     mkdir("prot_t5")
        #     wget.download('http://data.bioembeddings.com/public/embeddings/feature_models/t5/secstruct_checkpoint.pt',
        #                   "prot_t5/")

        # === lets leave the model in CPU to preserve GPU memory
        # self.device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
        self.device = 'cpu'
        print(f"Embedding running on {self.device}")

        if reference == 'prot_bert':
            self.latent_size = 1024
            self.model, self.tokenizer = self.get_bert_model(reference)
            self.bert = True
        else:
            self.latent_size = 1024
            self.model, self.tokenizer = self.get_T5_model(reference)
            self.bert = False

        # TODO: make sure that the meaning of '.' and '-' are correct
        # lets convert the alphabet to the latent space. We will need this to de-embed
        seq = {}
        self.lengths = {}
        seq['alphabet'] = ''.join(alphabet)
        # seq['reversed'] = tmp.join(alphabet)[::-1]

        self.latent_alphabet = np.zeros((len(alphabet), self.latent_size))
        for i, a in enumerate(alphabet):
            sequence = {
                'a1': f"_{a}_",
                'a2': f"X{a}X"
            }
            embedded = self.embed(sequence)
            # self.latent_alphabet[i, :] = embedded[1,1,:]
            self.latent_alphabet[i, :] = (embedded[0,1,:] + embedded[1,1,:])/2

        # embedded = self.embed(seq)
        # # average forward and reverse alphabets. Strip the very first and last tokens as it is the "start" and  "stop"
        # # self.latent_alphabet = (embedded[0, 1:-1, :] + np.flip(embedded[1, 1:-1, :], 0)) / 2
        # self.latent_alphabet = embedded[0, 1:-1, :]

    # Load encoder-part of ProtT5 in half-precision.
    # Load ProtT5 in half-precision (more specifically: the encoder-part of ProtT5-XL-U50)
    def get_T5_model(self, reference):
        model = T5EncoderModel.from_pretrained(f"Rostlab/{reference}")
        # === lets leave the model in CPU to preserve GPU memory
        # if torch.cuda.is_available():
        model = model.to(self.device)
        #     print("Transferred model to GPU")
        model = model.eval() # set model to evaluation model
        tokenizer = T5Tokenizer.from_pretrained(f'Rostlab/{reference}', do_lower_case=False)

        return model, tokenizer

    # bert model uses different class to load
    # does not produce a good result. Perhaps tokenizer is different?
    def get_bert_model(self, reference):
        model = BertModel.from_pretrained(f"Rostlab/{reference}")
        # === lets leave the model in CPU to preserve GPU memory
        # device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
        #     print("Transferred model to GPU")
        model = model.to(self.device)
        model = model.eval()         # set model to evaluation model
        tokenizer = BertTokenizer.from_pretrained(f'Rostlab/{reference}', do_lower_case=False)
        return model, tokenizer

    def embed(self, seqs):
        # memorise the lengths of chains
        self.max_len = 0
        for key, val in sorted(seqs.items()):
            val = str(val)
            # for rep in alphabet:
            #     val = str.replace(val, rep, '$')
            lv = len(val)
            # maximum seq lenght - will be the batch width
            self.lengths[key] = lv
            if lv > self.max_len:
                self.max_len = lv

        batch = list()
        for key, val in sorted(seqs.items()):
            # add spaces between letters (needed for tokenizer)
            batch.append(" ".join(val))

        # add_special_tokens adds extra token at the end of each sequence
        token_encoding = self.tokenizer.batch_encode_plus(batch, add_special_tokens=True, padding="longest", pad_to_max_length=True)
        input_ids = torch.tensor(token_encoding['input_ids']).to(self.device)
        attention_mask = torch.tensor(token_encoding['attention_mask']).to(self.device)

        try:
            with torch.no_grad():
                # returns: ( batch-size x max_seq_len_in_minibatch x embedding_dim )
                embedding_repr = self.model(input_ids, attention_mask=attention_mask)
        except RuntimeError:
            print("RuntimeError during embedding")

        # logits = np.zeros((len(seqs), self.max_len, self.latent_size))
        # for batch_idx, identifier in enumerate(seqs):  # for each protein in the current mini-batch
        #     s_len = self.lengths[identifier]
        #     # slice off padding --> batch-size x seq_len x embedding_dim
        #     emb = embedding_repr.last_hidden_state[batch_idx, :s_len].cpu().numpy()
        #     logits[batch_idx, 0:self.lengths[identifier], :] = emb


        # logits = embedding_repr.last_hidden_state[:, :self.max_len].cpu().numpy()
        logits = embedding_repr.last_hidden_state.cpu().numpy()

        # Remove padding ([PAD]) and special tokens ([CLS],[SEP]) that is added by Bert model
        # Bert and T5 have different tokenizers! Bert adds [PAD] at beginning and [CLS] at the end, while T5 only [CLS] at the end
        # features = []
        # for seq_num in range(len(logits)):
        #     # seq_len = (attention_mask[seq_num] == 1).sum()
        #     seq_len = self.max_len
        #     if self.bert:
        #         seq_emd = logits[seq_num][1:seq_len-1]
        #     else:
        #         seq_emd = logits[seq_num][0:seq_len-1]
        #     features.append(seq_emd)
        # we want to return padded array with the same length for each seq. So, we use maxlen instead of real len
        features = logits[:,]
        if self.bert:
            features = logits[:,1:-1,:]
        else:
            features = logits[:,0:-1,:]

        return features

    def de_embed(self, eee):
        seqs = {}

        # preserving the order!
        for k, chain in zip([0,1], ['H','L']):
            tmp = []

            for i in range(self.lengths[chain]):
                dists = []
                for a in self.latent_alphabet:
                    # dists.append(np.linalg.norm(a-eee[k, i, :]))   # euclidean distance - no good!
                    dists.append(distance.cosine(a, eee[k, i, :]))   # cosine distance

                match = np.argmin(dists)
                # print(f"{self.alphabet.unique_no_split_tokens[match]}", end='')
                tmp.append(alphabet[match])
            blank = ''
            seqs[chain] = blank.join(tmp)
            # print('')

        return seqs


if __name__ == '__main__':
    seqs = {
        'H': Seq('QVQLVESGGGLIQPGGSLRLSCAASGFIVSRNYMIWVRQAPGKGLEWVSVIYSGGSTFYA'),
        'L': Seq('EIVLTQSPGTLSLSPGERATLSCRASQSISSSYLAWYQQKPGQAPRLLIYGATSRATGTP'),
        # 'L': Seq(''.join(alphabet)),
    }

    # see https://huggingface.co/models?sort=downloads&search=rostlab
    # em = EmbedT5('prot_t5_xl_half_uniref50-enc')
    em = EmbedT5('prot_t5_xl_bfd')
    # em = EmbedT5('prot_bert')

    eee = em.embed(seqs)

    de_eee = em.de_embed(eee)
    # print([str(val) for key, val in seqs.items()])
    # print([val for key, val in de_eee.items()])
    sq1h = seqs['H']  # + ' ' + start_seq_spacers['L']
    sq2h = de_eee['H']  # + ' ' + np2seq_show(es.result.xfavorite)[1]
    sq1l = seqs['L']
    sq2l = de_eee['L']
    print(sq1h + " " + sq1l)
    # print(sq2h + " " + sq2l)
    highligted_h, changes_h = highlight_differences(sq1h, sq2h)
    highligted_l, changes_l = highlight_differences(sq1l, sq2l)
    print(f"{highligted_h} {highligted_l} Diff: {changes_h}+{changes_l}")

    # aaa =  dct(eee, type=2, n=800)
    # bbb = idct(aaa, type=2, n=1024)
    eee[:,:,800:] = 0
    de_eee = em.de_embed(eee)
    sq1h = seqs['H']  # + ' ' + start_seq_spacers['L']
    sq2h = de_eee['H']  # + ' ' + np2seq_show(es.result.xfavorite)[1]
    sq1l = seqs['L']
    sq2l = de_eee['L']
    print(sq1h + " " + sq1l)
    # print(sq2h + " " + sq2l)
    highligted_h, changes_h = highlight_differences(sq1h, sq2h)
    highligted_l, changes_l = highlight_differences(sq1l, sq2l)
    print(f"{sq2h} {sq2l} Diff: {changes_h}+{changes_l}")

    aaa =  dct(eee, type=3, n=800)
    bbb = idct(aaa, type=3, n=1024)
    de_eee = em.de_embed(bbb)
    sq1h = seqs['H']  # + ' ' + start_seq_spacers['L']
    sq2h = de_eee['H']  # + ' ' + np2seq_show(es.result.xfavorite)[1]
    sq1l = seqs['L']
    sq2l = de_eee['L']
    print(sq1h + " " + sq1l)
    # print(sq2h + " " + sq2l)
    highligted_h, changes_h = highlight_differences(sq1h, sq2h)
    highligted_l, changes_l = highlight_differences(sq1l, sq2l)
    print(f"{sq2h} {sq2l} Diff: {changes_h}+{changes_l}")
