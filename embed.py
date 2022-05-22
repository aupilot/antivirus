import numpy as np
from Bio.Seq import Seq
import torch
# pip install fair-esm
from esm import Alphabet, FastaBatchedDataset, ProteinBertModel, pretrained, MSATransformer

# inspired by examples from
# https://github.com/facebookresearch/esm

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

        # TODO: make sure that the meaning of '.' and '-' are correct
        # lets convert the alphabet to the latent space. We will need this to de-embed
        seq = {}
        tmp = ''
        seq['alphabet'] = tmp.join(self.alphabet.unique_no_split_tokens)
        # seq['reversed'] = tmp.join(self.alphabet.unique_no_split_tokens)[::-1]
        embedded = self.embed(seq)
        # average forward and reverse alphabets. Strip the very first and last tokens as it is the "start" and  "stop"
        # self.latent_alphabet = (embedded[0, 1:-1, :] + np.flip(embedded[1, 1:-1, :], 0)) / 2
        self.latent_alphabet = embedded[0, 1:-1, :]


    # TODO: this will encode 2 chains. In future lets encode all candidates in one batch
    def embed(self, seqs):
        dataset = SeqDataset.from_seqs(seqs)
        batches = dataset.get_batch_indices(4096, extra_toks_per_seq=1)
        data_loader = torch.utils.data.DataLoader(
            dataset, collate_fn=self.alphabet.get_batch_converter(), batch_sampler=batches, shuffle=False
        )

        with torch.no_grad():
            for batch_idx, (labels, strs, toks) in enumerate(data_loader):
                print(
                    f"Processing {batch_idx + 1} of {len(batches)} batches ({toks.size(0)} sequences)"
                )

                out = self.model(toks, repr_layers=[0], return_contacts=False)
                logits = out["logits"].cpu().numpy()

        return logits

    def de_embed(self, eee):
        batch_size = 2
        eee_shaped = np.reshape(eee, (batch_size, -1, self.latent_size))

        seqs = {}
        # the ESM dataset sorts data by size. We replaced it with our method get_batch_indices that does not sort, thus preserving the order
        for k, chain in zip([0,1], ['H','L']):
            tmp = []

            # we skip the first letter ad it is always <pad> or <cls>
            for i in range(1, eee_shaped.shape[1]):
                dists = []
                for a in self.latent_alphabet:
                    dists.append(np.linalg.norm(a-eee_shaped[k,i,:]))
                match = np.argmin(dists)
                # print(f"{self.alphabet.unique_no_split_tokens[match]}", end='')
                if match < 2:  # <pad> == 1 indicated the end (and sometimes start). So break on <cls> or <pad>
                    break
                # skip all other system tokens
                if match < 4:
                    continue
                if match > 30:
                    continue
                # we don't allow extra residues such as X B U Z O. Replace with '-'
                if match >= 24 and match <= 28:
                    tmp.append(self.alphabet.unique_no_split_tokens[30])
                    continue
                tmp.append(self.alphabet.unique_no_split_tokens[match])
            blank = ''
            seqs[chain] = blank.join(tmp)
            # print('')

        return seqs


# override the dataset to be initialised from Biopython sequence dict
class SeqDataset(FastaBatchedDataset):

    @classmethod
    def from_seqs(cls, seqs):

        sequence_labels, sequence_strs = [], []

        # we wanted top keep the order
        for key, val in sorted(seqs.items()):
            sequence_labels.append(key)
            sequence_strs.append(val)
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
            'H': Seq('QLVLTQSPSASASLGASVKLTCTLSSGHSNYA'),
            'L': Seq('GSSSGAERY'),
    }
    # dataset = SeqDataset.from_seqs(seqs)
    # print(len(dataset))

    em = Embed("esm1b_t33_650M_UR50S")
    eee = em.embed(seqs)

    de_eee = em.de_embed(eee.flatten())
    print(de_eee)



