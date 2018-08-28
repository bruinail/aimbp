function [note] = get_note(admission, index, data_dir)
A = load(sprintf('%s/adm_%d.mat',data_dir,admission));
note = A.data_notes.text{index};