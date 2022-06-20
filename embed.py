import numpy as np
from Bio.Seq import Seq
import torch
# pip install fair-esm
from esm import Alphabet, FastaBatchedDataset, ProteinBertModel, pretrained, MSATransformer
from scipy.spatial import distance

# inspired by examples from
# https://github.com/facebookresearch/esm
from helper import highlight_differences


class Embed(object):
    def __init__(self, model_file):
        self.model, self.alphabet = pretrained.load_model_and_alphabet(model_file)
        self.model.eval()
        if isinstance(self.model, MSATransformer):
            raise ValueError(
                "This script currently does not handle models with MSA input (MSA Transformer)."
            )

        self.latent_size = self.model.num_layers
        # === lets leave the model in CPU to preserve GPU memory
        # if torch.cuda.is_available():
        #     self.model = self.model.cuda()
        #     print("Transferred model to GPU")
        print("Embedding running on CPU")

        # TODO: make sure that the meaning of '.' and '-' are correct
        # lets convert the alphabet to the latent space. We will need this to de-embed
        seq = {}
        self.lengths = {}
        tmp = ''

        # self.latent_alphabet = np.zeros((len(self.alphabet.unique_no_split_tokens), self.latent_size))
        # for i, a in enumerate(self.alphabet.unique_no_split_tokens):
        #     sequence = {
        #         'a1': f"-{a}-",
        #         'a2': f"X{a}X"
        #     }
        #     embedded = self.embed(sequence)
        #     # self.latent_alphabet[i, :] = embedded[1,1,:]
        #     self.latent_alphabet[i, :] = (embedded[0,1,:] + embedded[1,1,:])/2

        seq['alphabet'] = tmp.join(self.alphabet.unique_no_split_tokens)
        # seq['reversed'] = tmp.join(self.alphabet.unique_no_split_tokens)[::-1]
        embedded = self.embed(seq)
        # average forward and reverse alphabets. Strip the very first and last tokens as it is the "start" and  "stop"
        # self.latent_alphabet = (embedded[0, 1:-1, :] + np.flip(embedded[1, 1:-1, :], 0)) / 2
        self.latent_alphabet = embedded[0, 1:-1, :]


    # TODO: this will encode 2 chains. In future lets encode all candidates in one batch
    def embed(self, seqs):
        # memorise the lengths of chains
        for key, val in sorted(seqs.items()):
            val = str(val)
            for rep in self.alphabet.unique_no_split_tokens:
                val = str.replace(val, rep, '$')
            self.lengths[key] = len(val)

        dataset = SeqDataset.from_seqs(seqs)
        batches = dataset.get_batch_indices(4096, extra_toks_per_seq=1)
        data_loader = torch.utils.data.DataLoader(
            dataset, collate_fn=self.alphabet.get_batch_converter(), batch_sampler=batches, shuffle=False
        )

        # attension: this only works for one batch with several sequences
        with torch.no_grad():
            for batch_idx, (labels, strs, toks) in enumerate(data_loader):
                print(f"Processing {batch_idx + 1} of {len(batches)} batches ({toks.size(0)} sequences)")
                out = self.model(toks, repr_layers=[], return_contacts=False)
                logits = out["logits"].cpu().numpy()

        return logits

    def de_embed(self, eee):
        # batch_size = 2
        # eee_shaped = np.reshape(eee, (batch_size, -1, self.latent_size))

        seqs = {}
        # the ESM dataset sorts data by size. We replaced it with our method get_batch_indices that does not sort, thus preserving the order
        for k, chain in zip([0,1], ['H','L']):
            tmp = []

            # we skip the first letter as it is always <pad> or <cls>. We also use memorised lengths as otherwise it will be aligned to max chain length
            for i in range(1, self.lengths[chain]+1):
            # for i in range(1, eee_shaped.shape[1]):
                dists = []
                for a in self.latent_alphabet:
                    # euq distance
                    # dists.append(np.linalg.norm(a-eee_shaped[k, i, :]))
                    # cosine distance
                    dists.append(distance.cosine(a, eee[k, i, :]))

                match = np.argmin(dists)
                # print(f"{self.alphabet.unique_no_split_tokens[match]}", end='')
                if match < 2:  # <pad> == 1 indicated the end (and sometimes start). So break on <cls> or <pad> # never got here
                    break
                # skip all other system tokens
                if match < 4:
                    continue
                if match > 30:
                    continue
                # if match == '<eos>':
                #     break  # never got here. Same with '.'
                # we don't allow extra residues such as X B U Z O. Also WTF is "."? Replace with all '-'
                if match >= 24 and match <= 29:
                    tmp.append('-')
                    continue
                tmp.append(self.alphabet.unique_no_split_tokens[match])
            blank = ''
            seqs[chain] = blank.join(tmp)
            # print('')

        return seqs
# esm alphabet
# ['<cls>', '<pad>', '<eos>', '<unk>', 'L', 'A', 'G', 'V', 'S', 'E', 'R', 'T', 'I', 'D', 'P', 'K', 'Q', 'N', 'F', 'Y', 'M', 'H', 'W', 'C', 'X', 'B', 'U', 'Z', 'O', '.', '-', '<null_1>', '<mask>']


# override the dataset to be initialised from Biopython sequence dict
class SeqDataset(FastaBatchedDataset):

    @classmethod
    def from_seqs(cls, seqs):

        sequence_labels, sequence_strs = [], []

        # we wanted top keep the order
        for key, val in sorted(seqs.items()):
            sequence_labels.append(key)
            # sequence_strs.append(val)
            seq_str = str(val)
            # for some reason the "-" or any other symbol does not encode spacers... too bad
            # if key == "H" or key == "L":
            #     seq_str = seq_str.replace("-", ".")
            sequence_strs.append(seq_str + ".")          # convert to strings because we have to use different dictionary - compatible with ESM
            # the below does not work. Insted we need to memorise the length of the chains in embed(). Not ideal...
            # sequence_strs.append(val + "<cls>")  # we adding the "eos" to help the embedding limit the sequence
        assert len(set(sequence_labels)) == len(
            sequence_labels
        ), "Found duplicate sequence labels"

        return cls(sequence_labels, sequence_strs)

    def get_batch_indices(self, toks_per_batch, extra_toks_per_seq=0):
        sizes = [(len(s), i) for i, s in enumerate(self.sequence_strs)]
        # sizes.sort()      # ===kir! preserve the order
        batches = []
        buf = []
        max_len = 0

        def _flush_current_buf():
            nonlocal max_len, buf
            if len(buf) == 0:
                return
            batches.append(buf)
            buf = []
            max_len = 0

        for sz, i in sizes:
            sz += extra_toks_per_seq
            if max(sz, max_len) * (len(buf) + 1) > toks_per_batch:
                _flush_current_buf()
            max_len = max(max_len, sz)
            buf.append(i)

        _flush_current_buf()
        return batches


if __name__ == '__main__':
    seqs = {
            'H': Seq('QVQLVESGGGLIQPGGSLRLSCAASGFIVSRNYMIWVRQAPGKGLEWVSVIYSGGSTFYA'),
            'L': Seq('EIVLTQSPGTLSLSPGERATLSCRASQSISSSYLAWYQQKPGQAPRLLIYGATSRATGTP'),
    }
    # dataset = SeqDataset.from_seqs(seqs)
    # print(len(dataset))

    em = Embed("esm1v_t33_650M_UR90S_5")
    eee = em.embed(seqs)
    de_eee = em.de_embed(eee)

    sq1h = seqs['H']  # + ' ' + start_seq_spacers['L']
    sq2h = de_eee['H']  # + ' ' + np2seq_show(es.result.xfavorite)[1]
    sq1l = seqs['L']
    sq2l = de_eee['L']
    print(sq1h + " " + sq1l)
    # print(sq2h + " " + sq2l)
    highligted_h, changes_h = highlight_differences(sq1h, sq2h)
    highligted_l, changes_l = highlight_differences(sq1l, sq2l)
    print(f"{highligted_h} {highligted_l} Diff: {changes_h}+{changes_l}")

    # print([str(val) for key, val in seqs.items()])
    # print([val for key, val in de_eee.items()])



