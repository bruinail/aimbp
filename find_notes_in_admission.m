function [found] = find_notes_in_admission(admission, search_regex, data_dir)
% Finds notes in a given admission that matches a search term

A = load(sprintf('%s/adm_%d.mat',data_dir,admission));
if(height(A.data_notes) == 0)
    found = cell(0,0);
    return;
end
found = cellfun(@(x) regexp(x,search_regex),A.data_notes.text,'UniformOutput',false);
