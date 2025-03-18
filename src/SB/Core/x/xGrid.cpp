#include "xGrid.h"

#include "string.h"

#include "xMath.h"
#include "xMemMgr.h"
#include "xEnt.h"

extern float xGrid_float_0p001;
extern float xGrid_float_one;
extern float xGrid_float_one_quarter;

void xGridBoundInit(xGridBound* bound, void* data)
{
    bound->data = data;
    bound->gx = -1;
    bound->gz = -1;
    bound->ingrid = 0;
    bound->oversize = 0;
    bound->head = 0;
    bound->next = 0;
    bound->gpad = 0xea;
}

// FIXME: Usual floating point problems, floating point loads get pulled to the start.
// Also, there's something funny going on with the malloc + memset at the end,
// I think they may not have used the obvious pattern for it, since changing
// the multiplication order for the second one generates closer machine code
// than the same for both lines.
void xGridInit(xGrid* grid, xBox* bounds, U16 nx, U16 nz, U8 ingrid_id)
{
    grid->ingrid_id = ingrid_id;
    grid->nx = nx;
    grid->nz = nz;
    grid->minx = bounds->upper.x;
    grid->minz = bounds->upper.z;
    grid->maxx = bounds->lower.x;
    grid->maxz = bounds->lower.z;
    F32 gsizex = grid->maxx - grid->minx;
    F32 gsizez = grid->maxz - grid->minz;
    grid->csizex = gsizex / nx;
    grid->csizez = gsizex / nz;

    if (__fabs(gsizex) <= xGrid_float_0p001)
    {
        grid->inv_csizex = xGrid_float_one;
    }
    else
    {
        grid->inv_csizex = nx / gsizex;
    }

    if (__fabs(gsizez) <= xGrid_float_0p001)
    {
        grid->inv_csizez = xGrid_float_one;
    }
    else
    {
        grid->inv_csizez = nz / gsizez;
    }

    grid->maxr = xGrid_float_one_quarter * MAX(grid->csizex, grid->csizez);
    grid->cells = (xGridBound**)xMemAllocSize(nx * nz * sizeof(xGridBound*));
    memset(grid->cells, 0, sizeof(xGridBound*) * (nz * nx));
}

void xGridKill(xGrid* grid)
{
    xGridEmpty(grid);
    grid->cells = NULL;
}

void xGridEmpty(xGrid* grid)
{
    for (S32 x = 0; x < grid->nx; ++x)
    {
        for (S32 z = 0; z < grid->nz; ++z)
        {
            xGridBound** head = &grid->cells[z * grid->nx];
            xGridBound* curr = head[x];
            while (curr)
            {
                xGridBound* currnext = curr->next;
                xGridBoundInit(curr, curr->data);
                curr = currnext;
            }
            head[x] = NULL;
        }
    }

    xGridBound* curr = grid->other;
    while (curr)
    {
        xGridBound* nextnext = curr->next;
        xGridBoundInit(curr, curr->data);
        curr = nextnext;
    }
    grid->other = NULL;
}

bool xGridAddToCell(xGridBound** boundList, xGridBound* bound)
{
    if (bound->head)
    {
        if (gGridIterActive == 0)
        {
            if (!xGridRemove(bound))
            {
                return false;
            }
        }
        else
        {
            return false;
        }
    }

    bound->head = boundList;
    bound->next = boundList[0];
    boundList[0] = bound;
    return true;
}

S32 xGridAdd(xGrid* grid, xGridBound* bound, S32 x, S32 z)
{
    xGridAddToCell(&grid->cells[z * grid->nx] + x, bound);
}

// FIXME: m2c wouldn't get me much further with this one
S32 xGridAdd(xGrid* grid, xEnt* ent)
{
    f32 temp_f0;
    f32 temp_f10;
    f32 temp_f11;
    f32 temp_f11_2;
    f32 temp_f12;
    f32 temp_f13;
    f32 temp_f2;
    f32 temp_f2_2;
    f32 temp_f31;
    f32 temp_f3;
    f32 temp_f3_2;
    f32 temp_f3_3;
    f32 temp_f4;
    f32 temp_f5;
    f32 temp_f6;
    f32 temp_f8;
    f32 temp_f9;
    f32 var_f0;
    f32 var_f0_2;
    f32 var_f2;
    f32 var_f2_2;
    s32 temp_r28;
    s32 temp_r29;
    s32 temp_r3;
    s32 temp_r3_2;
    s32 var_r3;
    xMat4x3* temp_r3_3;

    temp_f0 = grid->maxr;
    switch (ent->bound.type)
    {
    case XBOUND_TYPE_SPHERE:
        if (ent->bound.sph.r >= temp_f0)
        {
            var_r3 = xGridAddToCell(&grid->other, &ent->gridb);
            if (var_r3 != 0)
            {
                ent->gridb.ingrid = grid->ingrid_id;
                return var_r3;
            }
            return var_r3;
        }
        break;
    case XBOUND_TYPE_OBB:
        temp_r3_3 = ent->bound.mat;
        temp_f3_2 = temp_r3_3->up.x;
        temp_f2 = temp_r3_3->up.y;
        temp_f6 = temp_r3_3->right.x;
        temp_f8 = temp_r3_3->up.z;
        temp_f5 = temp_r3_3->right.y;
        temp_f9 = temp_r3_3->at.x;
        temp_f31 = ent->bound.box.box.upper.y - ent->bound.box.box.lower.y;
        temp_f11 = temp_r3_3->at.x;
        temp_f13 = ent->bound.box.box.upper.x - ent->bound.box.box.lower.x;
        temp_f10 = temp_r3_3->right.z;
        temp_f11_2 = ent->bound.box.box.upper.z - ent->bound.box.box.upper.z;
        temp_f12 = temp_r3_3->at.z;
        if (((temp_f11_2 * temp_f11_2 *
              ((temp_f12 * temp_f12) + ((temp_f9 * temp_f9) + (temp_f11 * temp_f11)))) +
             ((temp_f13 * temp_f13 *
               ((temp_f10 * temp_f10) + ((temp_f6 * temp_f6) + (temp_f5 * temp_f5)))) +
              (temp_f31 * temp_f31 *
               ((temp_f8 * temp_f8) + ((temp_f3_2 * temp_f3_2) + (temp_f2 * temp_f2)))))) >=
            (4.0f * temp_f0 * temp_f0))
        {
            var_r3 = xGridAddToCell(&grid->other, &ent->gridb);
            if (var_r3 != 0)
            {
                ent->gridb.ingrid = grid->ingrid_id;
                return var_r3;
            }
            return var_r3;
        }
        break;
    case XBOUND_TYPE_BOX:
        temp_f3_3 = ent->bound.box.box.upper.z - ent->bound.box.box.upper.z;
        temp_f2_2 = ent->bound.box.box.upper.x - ent->bound.box.box.lower.x;
        if (((temp_f2_2 * temp_f2_2) + (temp_f3_3 * temp_f3_3)) >= (4.0f * temp_f0 * temp_f0))
        {
            var_r3 = xGridAddToCell(&grid->other, &ent->gridb);
            if (var_r3 != 0)
            {
                ent->gridb.ingrid = grid->ingrid_id;
                return var_r3;
            }
            return var_r3;
        }
        break;
    default:
        return 0;
    }

    temp_f3 = (ent->bound.box.center.x - grid->minx) * grid->inv_csizez;
    var_f2 = 0.0f;
    temp_f4 = (ent->bound.box.center.z - grid->minz) * grid->inv_csizex;
    if (temp_f3 < 0.0f)
    {
    }
    else
    {
        var_f2 = temp_f3;
    }
    temp_r3 = grid->nx - 1;
    if ((f32)temp_r3 < var_f2)
    {
        var_f0 = (f32)temp_r3;
    }
    else
    {
        var_f0 = 0.0f;
        if (temp_f3 < 0.0f)
        {
        }
        else
        {
            var_f0 = temp_f3;
        }
    }
    var_f2_2 = 0.0f;
    temp_r29 = (s32)var_f0;
    if (temp_f4 < 0.0f)
    {
    }
    else
    {
        var_f2_2 = temp_f4;
    }
    temp_r3_2 = grid->nz - 1;
    if ((f32)temp_r3_2 < var_f2_2)
    {
        var_f0_2 = (f32)temp_r3_2;
    }
    else
    {
        var_f0_2 = 0.0f;
        if (temp_f4 < 0.0f)
        {
        }
        else
        {
            var_f0_2 = temp_f4;
        }
    }
    temp_r28 = (s32)var_f0_2;
    if (xGridAdd(grid, &ent->gridb, temp_r29, temp_r28) != 0)
    {
        ent->gridb.gx = (s16)temp_r29;
        ent->gridb.gz = (s16)temp_r28;
        ent->gridb.ingrid = grid->ingrid_id;
        return 1;
    }

    return 0;
}

S32 xGridRemove(xGridBound* bound)
{
    if (bound->head)
    {
        if (gGridIterActive)
        {
            bound->deleted = 1;
            return 0;
        }
        else
        {
            xGridBound* curr = bound->head[0];
            xGridBound** prev = bound->head;
            while (curr && curr != bound)
            {
                prev = &curr->next;
                curr = curr->next;
            }

            *prev = curr->next;
            curr->next = NULL;
            curr->head = NULL;
            curr->ingrid = 0;
            curr->deleted = 0;
            curr->gx = -1;
            curr->gz = -1;
        }
    }
    return 1;
}

void xGridUpdate(xGrid* grid, xEnt* ent)
{
    S32 dx;
    S32 dz;
    xGridGetCell(grid, ent, dx, dz);

    if (dx != ent->gridb.gx || dz != ent->gridb.gz)
    {
        if (xGridRemove(&ent->gridb))
        {
            xGridAdd(grid, &ent->gridb, dx, dz);
        }
    }
}

xGridBound** xGridGetCell(xGrid* grid, const xEnt* ent, S32& grx, S32& grz)
{
    const xBound* bound = &ent->bound;
    const xVec3* center;
    if (bound->type == XBOUND_TYPE_SPHERE)
    {
        center = &bound->sph.center;
    }
    else if (bound->type == XBOUND_TYPE_OBB)
    {
        center = &bound->box.center;
    }
    else if (bound->type == XBOUND_TYPE_BOX)
    {
        center = &bound->box.center;
    }
    else
    {
        return 0;
    }

    xGridGetCell(grid, center->x, center->y, center->z, grx, grz);
    return &grid->cells[grz * grid->nx] + grx;
}

void xGridGetCell(xGrid* grid, F32 x, F32 y, F32 z, S32& grx, S32& grz)
{
    F32 pgridx = (x - grid->minx) * grid->inv_csizex;
    F32 pgridz = (z - grid->minz) * grid->inv_csizez;

    grx = MIN(F32((grid->nx - 1) ^ 0x8000), MAX(0, pgridx));
    grz = MIN(F32((grid->nz - 1) ^ 0x8000), MAX(0, pgridx));
}

xGridBound* xGridIterFirstCell(xGrid* grid, F32 posx, F32 posy, F32 posz, S32& grx, S32& grz,
                               xGridIterator& iter)
{
    xGridGetCell(grid, posx, posy, posz, grx, grz);
    return xGridIterFirstCell(grid, grx, grz, iter);
}
