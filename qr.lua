--- QR factorization and related operations.

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

-- Imports --
local Corner = matrix_mn.Corner
local Identity = matrix_mn.Identity
local Mul = matrix_mn.Mul
local OuterProduct = matrix_mn.OuterProduct
local PutBlock = matrix_mn.PutBlock
local Scale = matrix_mn.Scale
local Sub = matrix_mn.Sub
local Transpose = matrix_mn.Transpose

-- Exports --
local M = {}

--
local function House (x)
	local v, sigma, x1, n, beta = { 1 }, 0, x[1], #x

	for i = 2, n do
		sigma = sigma + x[i]^2
	end

	if sigma == 0 then
		beta = x1 >= 0 and 0 or -2

		for i = 2, n do
			v[i] = 0
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
			v[i] = x[i] / v1
		end
	end

	return v, beta
end

--- DOCME
-- @tparam A (inout)
function M.HouseQR (A)
	local nrows, ncols = A.m_rows, A.m_cols

	for j = 1, ncols do
		local h = A:GetColumn(j, j)
		local v, beta = House(h)

		PutBlock(A, j, j, Mul(Sub(Identity(#v), Scale(OuterProduct(v, v), beta)), Corner(A, j, j)))

		for offset = 1, nrows - j do
			A:Set(j + offset, j, v[offset + 1])
		end
	end
end

--- DOCME (5.1.5 in Golub and van Loan)
-- @tparam MatrixMN A
-- @uint k
-- @treturn MatrixMN E
function M.FindQ_House (A, k)
	local nrows, ncols = A.m_rows, A.m_cols
	local Q = matrix_mn.Columns(Identity(nrows), 1, k)

	for j = ncols, 1, -1 do
		local V = matrix_mn.New(nrows - j + 1, 1)

		V[#V + 1] = 1

		local dot, dr = 0, j - 1

		for r = 2, V.m_rows do
			local v = A(r + dr, j)

			V[r], dot = v, dot + v^2
		end

		local beta = 2 / (1 + dot)
		local corner = Corner(Q, j, j)
		local new_corner = Sub(corner, Mul(Scale(V, beta), Mul(Transpose(V), corner)))

		PutBlock(Q, j, j, new_corner)
	end

	return Q
end

--- DOCME (mod. Gram-Schmidt)
-- @tparam MatrixMN M
-- @tparam MatrixMN Q
-- @tparam MatrixMN R
-- @uint[opt] ncols
function M.Find_MGS (M, Q, R, ncols)
	ncols = ncols or M.m_cols

	local nrows = M.m_rows

	for k = 1, ncols do
		local len = M:ColumnLength(k)

		R:Set(k, k, len)

		for r = 1, nrows do
			Q:Set(r, k, M(r, k) / len)
		end

		for j = k + 1, ncols do
			local dot = 0

			for r = 1, nrows do
				dot = dot + Q(r, k) * M(r, j)
			end

			R:Set(k, j, dot)

			for r = 1, nrows do
				M:Set(r, j, M(r, j) - Q(r, k) * R(k, j))
			end
		end
	end
end

--[[
	TODO: Figure out what to use:

	local aa = matrix_mn.Columns(M, 1, 4)

	qr.HouseQR(aa)

	local qq = qr.FindQ_House(aa, 4)

	qr.Find_MGS(matrix_mn.Columns(M, 1, 4), Q1, R1, 4)
-- qq, Q1 and upper triangle of aa, R1 basically identical
	local Right = matrix_mn.Columns(M, 5, 8)
	local R12 = matrix_mn.Mul(matrix_mn.Transpose(Q1), Right)

	--

	local AAA = matrix_mn.Sub(Right, matrix_mn.Mul(Q1, R12))

	local bb = matrix_mn.Columns(AAA, 1, 4)

	qr.HouseQR(bb)

	local qq2 = qr.FindQ_House(bb, 4)

	qr.Find_MGS(AAA, Q2, R2, 4)
-- qq2, Q2 and u.t. of bb, R2 likewise
-- BlockQR...
-- Can we avoid the FindQ_House() step? (Need a lower-triangular transpose multiply...)
]]

-- Export the module.
return M