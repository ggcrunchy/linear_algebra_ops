--- Operations dealing with Householder reflections.

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
local sqrt = math.sqrt

-- Modules --
local matrix_mn = require("numeric_types.matrix_mn")

-- Exports --
local M = {}

--- DOCME
function M.GetVector (x, out)
	--
	local sigma, n = 0, x:GetRowCount()

	if out then
		matrix_mn.Resize(out, 1, n)
	else
		out = matrix_mn.New(1, n)
	end

	for i = 2, n do
		sigma = sigma + x(i, 1)^2
	end

	--
	local x1, beta = x(1, 1)

	out:Set(1, 1, 1)

	if sigma == 0 then
		beta = x1 >= 0 and 0 or -2

		for i = 2, n do
			out:Set(1, i, 0)
		end
	else
		local mu, v1 = sqrt(x1^2 + sigma)

		if x1 <= 0 then
			v1 = x1 - mu
		else
			v1 = -sigma / (x1 + mu)
		end

		local v1sq = v1^2

		beta = 2 * v1sq / (sigma + v1sq)

		for i = 2, n do
			out:Set(1, i, x(i, 1) / v1)
		end
	end

	return beta, out
end

-- Export the module.
return M