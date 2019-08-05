function varargout = daqfnAns(fn, varargin)
% function varargout = daqfn(fn, varargin)
% (c) 2010 Hendrik Bluhm.  Please see LICENSE and COPYRIGHT information in plssetup.m.


status = uint32(0);
%fprintf('Calling %s\n',fn);
[varargout{1:nargout}] = calllib('ATSApi', ['Alazar', fn], varargin{:});

end