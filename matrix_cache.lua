--- Caching utilities for m-by-n matrices used by linear algebra algorithms.

--
-- Permission is hereby granted, free of charge, to any person obtaining
-- a copy of this software and associated documentation files (the
-- "Software"), to deal in the Software without restriction, including
-- without limitation the rights to use, copy, modify, merge, publish,
-- distribute, sublicense, and/or sell copies of the Software, and to
-- permit persons to whom the Software is furnished to do so, subject to
-- the following conditions:
--
-- The above copyright notice and this permission notice shall be
-- included in all copies or substantial portions of the Software.
--
-- THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
-- EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
-- MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
-- IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
-- CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
-- TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
-- SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
--
-- [ MIT license: http://www.opensource.org/licenses/mit-license.php ]
--

-- Standard library imports --
local remove = table.remove

-- Modules --
local collect = require("tektite_core.array.collect")
local matrix_mn = require("numeric_types.matrix_mn")
local wipe = require("tektite_core.array.wipe")

-- Exports --
local M = {}

-- --
local Cache = {}

--- DOCME
-- @treturn Matrix_MN V
function M.GetMatrix ()
	return remove(Cache) or matrix_mn.New(1, 1)
end

--- DOCME
-- @uint n
-- @treturn Matrix_MN V
function M.GetMatrix_N (n)
	for i = #Cache + 1, n do
		Cache[i] = matrix_mn.New(1, 1)
	end

	return wipe.UnpackAndWipeRange(Cache, #Cache - n + 1)
end

--- DOCME
-- @param ...
function M.Recache (...)
	collect.AppendArgs(Cache, ...)
end

-- Export the module.
return M