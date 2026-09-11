function wrote = ea_write_provenance(fpath, series, url, overwrite)
% EA_WRITE_PROVENANCE  Record, next to an accepted raw data file, WHICH
% series it is and where it came from.
%
% assemble_dataset copies these lines into ds.meta.sources, so a saved
% dataset always names its own inputs.
%
% THE overwrite FLAG IS THE POINT OF THIS FUNCTION
% ------------------------------------------------
% fetch_outcome_data writes provenance in two situations, and they must
% behave differently:
%   * it DOWNLOADED the file -- it knows the exact URL, so it writes
%     (overwrite = true);
%   * it found a valid file already on disk -- all it can honestly say
%     is "not recorded by this run", and writing that over a real source
%     URL recorded by an earlier fetch would DESTROY the provenance the
%     dataset is supposed to carry.  So it passes overwrite = false and
%     the existing sidecar is left alone.
%
% INPUTS
% ------
%   fpath     : the raw data file; the sidecar is <fpath>.source.txt.
%   series    : human-readable series identifier.
%   url       : source URL, or a note when there is none.
%   overwrite : default true; false leaves an existing sidecar alone.
%
% OUTPUT
% ------
%   wrote : true if the sidecar was (re)written, false if it was left
%           alone or could not be written.
%
% NOTES
% -----
% Provenance is a convenience, never a reason to fail a fetch: an I/O
% error here returns false rather than raising.

if nargin < 4, overwrite = true; end
sidecar = [fpath '.source.txt'];
wrote = false;
if ~overwrite && exist(sidecar, 'file') == 2
    return
end
try
    fid = fopen(sidecar, 'w');
    if fid <= 0, return; end
    fprintf(fid, 'file        = %s\n', basename_(fpath));
    fprintf(fid, 'series      = %s\n', series);
    fprintf(fid, 'source_url  = %s\n', url);
    fprintf(fid, 'recorded_at = %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));  %#ok<TNOW1,DATST>
    fclose(fid);
    wrote = true;
catch
    % see NOTES
end
end

% -------------------------------------------------------------------------
function b = basename_(p)
[~, n, e] = fileparts(p);
b = [n e];
end
