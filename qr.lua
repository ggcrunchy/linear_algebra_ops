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
local floor = math.floor
local remove = table.remove
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

-- Cached module references --
local _BlockQR_
local _Multiply_TranposeHouseholder_
local _ThinQR_HouseholderColumns_

-- Exports --
local M = {}

-- --
local Cache = {}

--
local function GetMatrix ()
	return remove(Cache) or matrix_mn.New(1, 1)
end

--
local function Recache (a, b, c, d)
	Cache[#Cache + 1] = a
	Cache[#Cache + 1] = b
	Cache[#Cache + 1] = c
	Cache[#Cache + 1] = d
end

--
local function AuxBlockQR (A, col, n, nb, out)
	if n <= nb then
		_ThinQR_HouseholderColumns_(A, col, n, out)
	else
		local n1 = floor(.5 * n)

		AuxBlockQR(A, 1, n1, nb, out)
		-- [Q1, R11] = BlockQR(A(:, 1:n1), n1, nb)

--		_Multiply_TranposeHouseholder_(Q1, A, R12)
		-- R12 = Q1^T * A(:, n1 + 1:n)

		-- A(:, n1 + 1:n) = A(:, n1 + 1:n) - Q1 * R12
		-- [Q2, R22] = BlockQR(A(:, n1 + 1:n), n - n1, nb)
		-- Q = [Q1|Q2], R = [R11 R12]
		--					[  0 R22]
	end
end

--
local function AuxOut (A, out)
	if out then
		out:Resize(A:GetDims())

		return out
	else
		return A
	end
end

--- DOCME
-- @tparam MatrixMN A
-- @uint nb
-- @tparam[opt=A] MatrixMN out
function M.BlockQR (A, nb, out)
	AuxBlockQR(A, 1, A:GetColumnCount(), nb, AuxOut(A, out))
end

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
-- @tparam MatrixMN A (inout)
-- @tparam[opt=A] MatrixMN out
function M.HouseQR (A, out)
	_ThinQR_HouseholderColumns_(A, 1, A:GetColumnCount(), out)
end

--- DOCME (5.1.5 in Golub and van Loan)
-- @tparam MatrixMN A
-- @uint k
-- @treturn MatrixMN E
function M.FindQ_House (A, k) -- rename: HouseholderToQ, add *ToQR variant...
	local nrows, ncols = A:GetDims()
	local Q = matrix_mn.Columns(Identity(nrows), 1, k)

	for j = ncols, 1, -1 do
		local vrows = nrows - j + 1
		local V = matrix_mn.New(vrows, 1)

		V[#V + 1] = 1

		local dot, dr = 0, j - 1

		for r = 2, vrows do
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
-- @tparam MatrixMN A
-- @tparam MatrixMN Q
-- @tparam MatrixMN R
-- @uint[opt] ncols
function M.Find_MGS (A, Q, R, ncols)
	ncols = ncols or A:GetColumnCount()

	local nrows = A:GetRowCount()

	for k = 1, ncols do
		local len = A:ColumnLength(k)

		R:Set(k, k, len)

		for r = 1, nrows do
			Q:Set(r, k, A(r, k) / len)
		end

		for j = k + 1, ncols do
			local dot = 0

			for r = 1, nrows do
				dot = dot + Q(r, k) * A(r, j)
			end

			R:Set(k, j, dot)

			for r = 1, nrows do
				A:Update(r, j, -Q(r, k) * R(k, j))
			end
		end
	end
end

-- --
local ColumnOpts = { to = 2 }

--- DOCME
function M.Multiply_TranposeHouseholder (H, C, out)
	out = AuxOut(C, out)

	--
	local corner, vc, v, op = GetMatrix(), GetMatrix(), GetMatrix(), GetMatrix()

	ColumnOpts.out = v

	for j = 1, C:GetColumnCount() do
		--
		ColumnOpts.from = j + 1

		H:GetColumn(j, ColumnOpts)
		v:Set(1, 1, 1)

		--
		Corner(C, j, 1, corner)
		Transpose(v, v)
		Mul(v, corner, vc)

		--
		Transpose(v, v)
		Scale(v, 2 / matrix_mn.FrobeniusNormSquared(v), v)
		PutBlock(out, j, 1, Sub(corner, OuterProduct(v, vc, op), corner))

		--
		C = out
	end

	Recache(corner, vc, v, op)

	ColumnOpts.out = nil
end
-- ^^^ TODO: Columns variant...

--- DOCME
-- @tparam MatrixMN A
-- @tparam Vector X
-- @tparam Vector B
-- @uint nb
-- @tparam[opt=A] MatrixMN out
function M.Solve_Householder (A, X, B, nb, out)
	_BlockQR_(A, nb, out)

	-- Back-sub!
	-- Q^t*XX
end

--- DOCME
-- @tparam MatrixMN A (inout)
-- @uint col
-- @uint n
-- @tparam[opt=A] MatrixMN out
function M.ThinQR_HouseholderColumns (A, col, n, out)
	out = AuxOut(A, out)

	--
	local dc, nrows = col - 1, A:GetRowCount()

	for j = 1, n do
		local ci = j + dc
		local h = A:GetColumn(ci, j)
		local v, beta = House(h)

		PutBlock(out, ci, j, Mul(Sub(Identity(#v), Scale(OuterProduct(v, v), beta)), Corner(A, ci, j)))

		for offset = 1, nrows - j do
			out:Set(j + offset, ci, v[offset + 1])
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

-- Cache module members.
_BlockQR_ = M.BlockQR
_Multiply_TranposeHouseholder_ = M.Multiply_TranposeHouseholder
_ThinQR_HouseholderColumns_ = M.ThinQR_HouseholderColumns

-- Export the module.
return M